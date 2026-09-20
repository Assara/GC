#pragma once
#include "block_minimal_generator.hpp"
#include <bit>

namespace VectorSpace::block_wiedemann_detail {
// PM-Basis: solve the first half, form the shifted residual, solve the
// second half with the updated shift, and multiply the two bases.
// Polynomial matrices are coefficient-major contiguous slabs. All scratch
// comes from arenas allocated at construction; recursion never allocates.
template<Field K>
class divide_conquer_generator {
    struct arena {
        OwnedArray<K> data;
        std::size_t used=0;
        explicit arena(std::size_t n):data(n) {}
        K* take(std::size_t n) {
            if(n>data.size()-used)throw std::length_error("PM-Basis arena exhausted");
            auto* p=data.data()+used;used+=n;std::fill_n(p,n,K{});return p;
        }
    };
    const packed_moments<K>& moments_;
    std::size_t b_,m_,training_,leaf_,multiply_leaf_;
    int threads_;
    minimal_generator_state<K> validation_;
    OwnedArray<K> input_,basis_;
    OwnedArray<std::size_t> indices_;
    arena scratch_;
    enum slot : std::size_t {order_s, input_s, output_s, shift_s, degrees_s, indices_s,
        mark_s, left_s, right_s, residual_s, phase_s, step_s, slots};
    OwnedArray<std::size_t> frames_;
    std::size_t depth_=0,completed_=0;
    bool initialized_=false,finished_=false;
    using Progress=typename minimal_generator_state<K>::Progress;
    static std::size_t add(std::size_t a,std::size_t b) {
        if(b>std::numeric_limits<std::size_t>::max()-a)throw std::length_error("PM-Basis size overflow");
        return a+b;
    }
    static std::size_t mul(std::size_t a,std::size_t b){return recurrence_product(a,b);}
    static std::size_t padded(std::size_t n) {
        if(n>std::bit_floor(std::numeric_limits<std::size_t>::max()))throw std::length_error("polynomial size overflow");
        return std::bit_ceil(n);
    }
    static std::size_t karatsuba_space(std::size_t n,std::size_t a,std::size_t b,std::size_t c,std::size_t cutoff) {
        if(n<=cutoff)return 0;
        return add(mul(n/2,add(a,b)),add(mul(n,c),karatsuba_space(n/2,a,b,c,cutoff)));
    }
    static std::size_t product_space(std::size_t x,std::size_t y,std::size_t r,std::size_t k,std::size_t c,std::size_t cutoff) {
        const auto n=padded(std::max(x,y)),a=mul(r,k),b=mul(k,c),z=mul(r,c);
        return add(mul(n,add(add(a,b),mul(2,z))),karatsuba_space(n,a,b,z,cutoff));
    }
    static std::size_t workspace(std::size_t order,std::size_t m,std::size_t b,std::size_t leaf,std::size_t cutoff) {
        if(order<=leaf)return mul(m,b);
        const auto l=order/2,r=order-l;
        const auto left=mul(l+1,mul(m,m)),right=mul(r+1,mul(m,m)),residual=mul(r,mul(m,b));
        const auto locals=add(add(left,right),residual);
        return add(locals,std::max({workspace(l,m,b,leaf,cutoff),workspace(r,m,b,leaf,cutoff),
            product_space(l+1,order,m,m,b,cutoff),product_space(r+1,l+1,m,m,m,cutoff)}));
    }
    // Multiplication over the matrix coefficient ring; operand order is preserved.
    void karatsuba(const K* a,const K* b,K* out,std::size_t n,std::size_t r,std::size_t k,std::size_t c) {
        const auto ac=r*k,bc=k*c,cc=r*c;
        std::fill_n(out,2*n*cc,K{});
        if(n<=multiply_leaf_) {
            const bool parallel=threads_>1 && n*n*r*k*c>=32768;
            (void)parallel;
#pragma omp parallel for schedule(static) num_threads(threads_) if(parallel)
            for(std::size_t i=0;i<r;++i)
                for(std::size_t x=0;x<n;++x)for(std::size_t y=0;y<n;++y)
                    for(std::size_t j=0;j<k;++j) {
                        const auto scalar=a[x*ac+i*k+j];if(scalar==K{})continue;
                        for(std::size_t z=0;z<c;++z)out[(x+y)*cc+i*c+z]+=scalar*b[y*bc+j*c+z];
                    }
            return;
        }
        const auto h=n/2,mark=scratch_.used;
        auto* sa=scratch_.take(h*ac);auto* sb=scratch_.take(h*bc);auto* cross=scratch_.take(n*cc);
        for(std::size_t i=0;i<h*ac;++i)sa[i]=a[i]+a[h*ac+i];
        for(std::size_t i=0;i<h*bc;++i)sb[i]=b[i]+b[h*bc+i];
        karatsuba(a,b,out,h,r,k,c);
        karatsuba(a+h*ac,b+h*bc,out+n*cc,h,r,k,c);
        karatsuba(sa,sb,cross,h,r,k,c);
        for(std::size_t i=0;i<n*cc;++i)cross[i]-=out[i]+out[n*cc+i];
        for(std::size_t i=0;i<n*cc;++i)out[h*cc+i]+=cross[i];
        scratch_.used=mark;
    }
    // Copy only the requested coefficient interval (middle product for residuals).
    void product(const K* a,std::size_t na,const K* b,std::size_t nb,K* out,
        std::size_t first,std::size_t count,std::size_t r,std::size_t k,std::size_t c) {
        const auto mark=scratch_.used,n=padded(std::max(na,nb));
        auto* aa=scratch_.take(n*r*k);auto* bb=scratch_.take(n*k*c);auto* result=scratch_.take(2*n*r*c);
        std::copy_n(a,na*r*k,aa);std::copy_n(b,nb*k*c,bb);
        karatsuba(aa,bb,result,n,r,k,c);
        std::fill_n(out,count*r*c,K{});
        if(first<2*n)std::copy_n(result+first*r*c,std::min(count,2*n-first)*r*c,out);
        scratch_.used=mark;
    }
    static std::size_t polynomial_length(const K* p,std::size_t length,std::size_t block) {
        while(length>1 && std::all_of(p+(length-1)*block,p+length*block,[](K x){return x==K{};}))--length;
        return length;
    }
    std::size_t arena_offset() const {return input_.size()+basis_.size();}
    K* coefficients(std::size_t offset) {
        if(offset<input_.size())return input_.data()+offset;
        if(offset<arena_offset())return basis_.data()+offset-input_.size();
        return scratch_.data.data()+offset-arena_offset();
    }
    std::size_t take_offset(std::size_t size) {
        const auto offset=arena_offset()+scratch_.used;scratch_.take(size);return offset;
    }
    void push(std::size_t order,std::size_t input,std::size_t output,
        std::size_t shift,std::size_t degrees,std::size_t indices) {
        if((depth_+1)*slots>frames_.size())throw std::length_error("PM-Basis frame capacity exhausted");
        auto* frame=frames_.data()+depth_++*slots;std::fill_n(frame,slots,0);
        frame[order_s]=order;frame[input_s]=input;frame[output_s]=output;
        frame[shift_s]=shift;frame[degrees_s]=degrees;frame[indices_s]=indices;
        frame[mark_s]=scratch_.used;
    }
    void prepare_input() {
        std::fill(input_.begin(),input_.end(),K{});
        for(std::size_t t=0;t<training_;++t)std::copy_n(moments_[t].data(),b_*b_,input_.data()+t*m_*b_);
        for(std::size_t i=0;i<b_;++i)input_[b_*b_+i*b_+i]=1;
    }
    void initialize() {
        prepare_input();
        for(std::size_t i=0;i<m_;++i)indices_[i]=i>=b_;
        push(training_,0,input_.size(),0,m_,2*m_);initialized_=true;
    }
    void leaf_step(std::size_t* frame) {
        const auto t=frame[step_s],mm=m_*m_,mb=m_*b_,mark=scratch_.used;
        auto* out=coefficients(frame[output_s]);const auto* f=coefficients(frame[input_s]);
        auto* degrees=indices_.data()+frame[degrees_s];
        auto* permutation=indices_.data()+frame[indices_s];auto* pivots=permutation+m_;auto* columns=pivots+m_;
        auto* discrepancy=scratch_.take(mb);
                std::fill_n(discrepancy,mb,K{});
                for(std::size_t d=0;d<=t;++d)for(std::size_t r=0;r<m_;++r)
                    for(std::size_t k=0;k<m_;++k) {
                        const auto scalar=out[d*mm+r*m_+k];if(scalar==K{})continue;
                        for(std::size_t j=0;j<b_;++j)discrepancy[r*b_+j]+=scalar*f[(t-d)*mb+k*b_+j];
                    }
                std::iota(permutation,permutation+m_,0);
                std::sort(permutation,permutation+m_,[&](auto a,auto b){return degrees[a]!=degrees[b]?degrees[a]<degrees[b]:a<b;});
                std::size_t count=0;
                for(std::size_t j=0;j<m_;++j) {
                    const auto r=permutation[j];
                    for(std::size_t i=0;i<count;++i) {
                        const auto p=pivots[i],c=columns[i];
                        if(discrepancy[r*b_+c]==K{})continue;
                        const auto factor=discrepancy[r*b_+c]*discrepancy[p*b_+c].inv();
                        for(std::size_t z=0;z<b_;++z)discrepancy[r*b_+z]-=factor*discrepancy[p*b_+z];
                        for(std::size_t d=0;d<=t;++d)for(std::size_t z=0;z<m_;++z)
                            out[d*mm+r*m_+z]-=factor*out[d*mm+p*m_+z];
                    }
                    std::size_t c=0;while(c<b_ && discrepancy[r*b_+c]==K{})++c;
                    if(c<b_){pivots[count]=r;columns[count++]=c;}
                }
                for(std::size_t i=0;i<count;++i) {
                    const auto r=pivots[i];
                    for(std::size_t d=t+1;d>0;--d)std::copy_n(out+(d-1)*mm+r*m_,m_,out+d*mm+r*m_);
                    std::fill_n(out+r*m_,m_,K{});++degrees[r];
                }

        scratch_.used=mark;++frame[step_s];++completed_;
    }
    void import_basis() {
        // Extraction and validation are unchanged and restartable from this basis.
        validation_.processed=training_;
        const auto* degrees=indices_.data()+m_;
        for(std::size_t r=0;r<m_;++r) {
            std::size_t length=training_+1;
            while(length>1 && std::all_of(basis_.data()+(length-1)*m_*m_+r*m_,
                basis_.data()+(length-1)*m_*m_+(r+1)*m_,[](K x){return x==K{};}))--length;
            validation_.lengths[r]=length;validation_.degrees[r]=degrees[r];
            for(std::size_t d=0;d<length;++d)
                std::copy_n(basis_.data()+d*m_*m_+r*m_,m_,validation_.row_data(r)+(length-1-d)*m_);
        }
    }
    void validate_frames() const {
        std::size_t expected_order=training_,input=0,output=input_.size(),shift=0,degrees=m_,index=2*m_,mark=0;
        std::size_t counted=0;
        for(std::size_t d=0;d<depth_;++d) {
            const auto* f=frames_.data()+d*slots;
            if(f[order_s]!=expected_order || f[input_s]!=input || f[output_s]!=output
                || f[shift_s]!=shift || f[degrees_s]!=degrees || f[indices_s]!=index || f[mark_s]!=mark)
                throw std::runtime_error("invalid PM-Basis checkpoint frame offsets");
            const auto phase=f[phase_s],l=expected_order/2,r=expected_order-l;
            if(expected_order<=leaf_) {
                if(d+1!=depth_ || (phase!=0 && phase!=4) || f[step_s]>=expected_order
                    || (phase==0 && f[step_s]) || index+3*m_>indices_.size())
                    throw std::runtime_error("invalid PM-Basis leaf frame");
                counted+=f[step_s];
            } else if(phase==0) {
                if(d+1!=depth_ || f[step_s])throw std::runtime_error("invalid PM-Basis pending frame");
            } else {
                if(phase>3 || f[step_s] || index+m_>indices_.size())throw std::runtime_error("invalid PM-Basis internal frame");
                const auto left=arena_offset()+mark,right=left+(l+1)*m_*m_,residual=right+(r+1)*m_*m_;
                if(f[left_s]!=left || f[right_s]!=right || f[residual_s]!=residual)
                    throw std::runtime_error("invalid PM-Basis arena offsets");
                mark=residual+r*m_*b_-arena_offset();
                if(mark>scratch_.used)throw std::runtime_error("invalid PM-Basis live arena");
                if(d+1<depth_) {
                    if(phase==1){expected_order=l;output=left;degrees=index;}
                    else if(phase==3){counted+=l;expected_order=r;input=residual;output=right;shift=index;}
                    else throw std::runtime_error("invalid PM-Basis child phase");
                    index+=m_;
                } else counted+=phase==3 ? expected_order : l;
            }
        }
        if(depth_ && (mark!=scratch_.used || counted!=completed_))throw std::runtime_error("inconsistent PM-Basis progress");
    }
    std::uint64_t sequence_identity() const {
        // Identity guard, not a cryptographic operator/data certificate.
        std::uint64_t hash=14695981039346656037ULL;
        const auto* bytes=reinterpret_cast<const unsigned char*>(moments_.data());
        const auto count=mul(mul(moments_.size(),b_*b_),sizeof(K));
        for(std::size_t i=0;i<count;++i){hash^=bytes[i];hash*=1099511628211ULL;}
        return hash;
    }

public:
    divide_conquer_generator(const packed_moments<K>& moments,std::size_t training,
        int threads=1,std::size_t leaf=16,std::size_t multiply_leaf=8)
        :moments_(moments),b_(moments.block_size()),m_(mul(2,b_)),training_(training),leaf_(leaf),multiply_leaf_(multiply_leaf),threads_(threads),
         validation_(moments,b_,threads),input_(mul(training,mul(m_,b_))),basis_(mul(add(training,1),mul(m_,m_))),
         indices_(mul(add(std::bit_width(training),6),m_)),
         scratch_(leaf && multiply_leaf ? workspace(training,m_,b_,leaf,multiply_leaf) : 0),
         frames_(mul(add(std::bit_width(training),1),slots)) {
        if(!training || training>moments.size() || !leaf || !multiply_leaf || threads<1)
            throw std::invalid_argument("invalid PM-Basis configuration");
    }
    // One stable transition: frame setup, one leaf term, or one full product.
    // A callback may save and throw to pause; all offsets already describe the
    // next operation. Karatsuba's temporary recursion is never serialized.
    bool step() {
        if(finished_)return false;
        if(!initialized_)initialize();
        auto* frame=frames_.data()+(depth_-1)*slots;
        const auto order=frame[order_s],l=order/2,r=order-l,mm=m_*m_,mb=m_*b_;
        auto* out=coefficients(frame[output_s]);
        switch(frame[phase_s]) {
        case 0:
            std::fill_n(out,(order+1)*mm,K{});
            if(order<=leaf_) {
                std::copy_n(indices_.data()+frame[shift_s],m_,indices_.data()+frame[degrees_s]);
                for(std::size_t i=0;i<m_;++i)out[i*m_+i]=1;
                frame[phase_s]=4;
            } else {
                frame[left_s]=take_offset((l+1)*mm);frame[right_s]=take_offset((r+1)*mm);
                frame[residual_s]=take_offset(r*mb);frame[phase_s]=1;
                push(l,frame[input_s],frame[left_s],frame[shift_s],frame[indices_s],frame[indices_s]+m_);
            }
            break;
        case 1: {
            const auto* left=coefficients(frame[left_s]);
            product(left,polynomial_length(left,l+1,mm),coefficients(frame[input_s]),order,
                coefficients(frame[residual_s]),l,r,m_,m_,b_);
            frame[phase_s]=2;break;
        }
        case 2:
            frame[phase_s]=3;
            push(r,frame[residual_s],frame[right_s],frame[indices_s],frame[degrees_s],frame[indices_s]+m_);
            break;
        case 3: {
            const auto* left=coefficients(frame[left_s]);const auto* right=coefficients(frame[right_s]);
            product(right,polynomial_length(right,r+1,mm),left,polynomial_length(left,l+1,mm),out,0,order+1,m_,m_,m_);
            scratch_.used=frame[mark_s];--depth_;break;
        }
        case 4:
            leaf_step(frame);
            if(frame[step_s]==order){scratch_.used=frame[mark_s];--depth_;}
            break;
        default:throw std::runtime_error("invalid PM-Basis frame phase");
        }
        if(!depth_){finished_=true;import_basis();}
        return !finished_;
    }
    void process(const Progress& progress={}) {
        while(!finished_) {step();if(progress)progress(completed_,training_);}
    }
    bool complete() const {return finished_;}
    std::size_t processed_terms() const {return completed_;}
    void save(serialization::archive_output& out) const {
        out.section("pm-basis-state",[&] {
            out.word(1);out.text(K::name());out.word(K::characteristic());out.word(sizeof(K));
            out.word(b_);out.word(training_);out.word(leaf_);out.word(multiply_leaf_);
            out.word(moments_.size());out.word(sequence_identity());out.word(initialized_);out.word(finished_);
            out.word(completed_);out.word(depth_);out.word(scratch_.used);
            if(initialized_) {
                out.integers<std::size_t>("work-stack",std::span<const std::size_t>(frames_.data(),depth_*slots));
                out.integers<std::size_t>("degree-workspace",indices_);
                out.block<K>("partial-basis",basis_,(training_+1)*m_,m_);
                out.block<K>("live-arena",std::span<const K>(scratch_.data.data(),scratch_.used),scratch_.used,1);
            }
        });
    }
    void load(serialization::archive_input& in) {
        in.section("pm-basis-state",[&] {
            in.expect(1);in.text(K::name());in.expect(K::characteristic());in.expect(sizeof(K));
            in.expect(b_);in.expect(training_);in.expect(leaf_);in.expect(multiply_leaf_);
            in.expect(moments_.size());in.expect(sequence_identity());
            const auto initialized=in.word(),finished=in.word();
            completed_=in.word();depth_=in.word();scratch_.used=in.word();
            if(initialized>1 || finished>1 || completed_>training_ || depth_>frames_.size()/slots
                || scratch_.used>scratch_.data.size() || (finished && (!initialized || depth_ || scratch_.used || completed_!=training_))
                || (!initialized && (finished || completed_ || depth_ || scratch_.used))
                || (initialized && !finished && !depth_))throw std::runtime_error("invalid PM-Basis checkpoint progress");
            initialized_=initialized;finished_=finished;
            if(initialized_) {
                in.integers<std::size_t>("work-stack",std::span<std::size_t>(frames_.data(),depth_*slots));
                in.integers<std::size_t>("degree-workspace",indices_);
                in.block<K>("partial-basis",basis_,(training_+1)*m_,m_);
                in.block<K>("live-arena",std::span<K>(scratch_.data.data(),scratch_.used),scratch_.used,1);
                validate_frames();prepare_input();
                if(finished_)import_basis();
            }
        });
    }
    bool generator(packed_generator<K>& result,const Progress& progress={}) {
        if(!finished_)throw std::logic_error("PM-Basis construction is incomplete");
        return validation_.generator(result,progress);
    }
    std::size_t allocated_bytes() const {
        return validation_.allocated_bytes()+(input_.size()+basis_.size()+scratch_.data.size())*sizeof(K)+(indices_.size()+frames_.size())*sizeof(std::size_t);
    }
};
}
