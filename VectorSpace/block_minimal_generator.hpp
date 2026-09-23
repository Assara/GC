#pragma once
#include <algorithm>
#include <numeric>
#include <functional>
#include <optional>
#include <vector>
#include "Field.hpp"
#include "packed_recurrence.hpp"

namespace VectorSpace::block_wiedemann_detail {
template <Field K>
std::size_t dense_rank(std::span<K> a, std::size_t n) {
    std::size_t rank = 0;
    for (std::size_t c = 0; c < n; ++c) {
        auto r = rank;
        while (r < n && a[r*n+c] == K{}) ++r;
        if (r == n) continue;
        for (std::size_t j = 0; j < n; ++j) std::swap(a[r*n+j], a[rank*n+j]);
        const K inverse = a[rank*n+c].inv();
        for (std::size_t i = rank+1; i < n; ++i) {
            const K factor = a[i*n+c] * inverse;
            for (std::size_t j = c; j < n; ++j) a[i*n+j] -= factor*a[rank*n+j];
        }
        ++rank;
    }
    return rank;
}

// Incremental shifted order basis of [S(z); I], shift (0,...,0,1,...,1).
// Like scalar Berlekamp--Massey, keep recurrence state when more moments arrive.
// Each polynomial row is flat and stored in descending powers: multiplying by
// z appends a zero block instead of moving all coefficients or allocating one
// vector per coefficient.
template <Field K>
class minimal_generator_state {
    template<Field> friend class divide_conquer_generator;
    const packed_moments<K>& moments;
    std::size_t b,width,processed=0;
    int threads;
    OwnedArray<K> basis,discrepancy,constant,leading;
    OwnedArray<std::size_t> lengths,degrees,order,selected,selected_degrees,validation_prefix;
    bool validation_active=false;
    struct pivot {std::size_t row,column;K inverse;};
    OwnedArray<pivot> pivots;
    std::size_t stride() const {return basis.size()/width;}
    K* row_data(std::size_t r) {return basis.data()+r*stride();}
    const K* row_data(std::size_t r) const {return basis.data()+r*stride();}
    std::size_t length(std::size_t r) const {return lengths[r];}
    const K* coefficient(std::size_t r,std::size_t d) const {
        return row_data(r)+(length(r)-1-d)*width;
    }
public:
    minimal_generator_state(const packed_moments<K>& sequence,std::size_t block,int workers=1)
        :moments(sequence),b(block),width(2*b),threads(workers),
         basis(recurrence_product(recurrence_product(width,width),sequence.capacity()+1)),
         discrepancy(width*b),constant(b*b),leading(b*b),lengths(width),degrees(width),order(width),
         selected(b),selected_degrees(b),validation_prefix(width),pivots(b) {
        if(!block || block!=sequence.block_size() || workers<1)throw std::invalid_argument("invalid recurrence dimensions");
        for(std::size_t r=0;r<width;++r){lengths[r]=1;row_data(r)[r]=1;degrees[r]=r>=b;}
    }
    std::size_t processed_terms() const {return processed;}
    std::size_t allocated_bytes() const {
        return (basis.size()+discrepancy.size()+constant.size()+leading.size())*sizeof(K)
            +(lengths.size()+degrees.size()+order.size()+selected.size()+selected_degrees.size()+validation_prefix.size())*sizeof(std::size_t)
            +pivots.size()*sizeof(pivot);
    }
    void save(serialization::archive_output& out) const {
        out.section("recurrence-state-v2",[&] {
            out.word(b);out.word(processed);
            out.word(validation_active);out.integers<std::size_t>("validation-prefix",validation_prefix);
            for(std::size_t r=0;r<width;++r) {
                out.word(length(r));out.word(degrees[r]);
                out.block<K>("polynomial-row",std::span<const K>(row_data(r),length(r)*width),length(r),width);
            }
        });
    }
    void load(serialization::archive_input& in) {
        in.section("recurrence-state-v2",[&] {
            in.expect(b);const auto terms=in.word();
            if(terms>moments.size())throw std::runtime_error("checkpoint recurrence exceeds available sequence");
            const auto active=in.word();if(active>1)throw std::runtime_error("invalid validation state");
            validation_active=active;
            in.integers<std::size_t>("validation-prefix",validation_prefix);
            for(auto t:validation_prefix)if(t>moments.size())throw std::runtime_error("invalid validation progress");
            for(std::size_t r=0;r<width;++r) {
                const auto len=in.word(),degree=in.word();
                if(!len || len>terms+1 || len>stride()/width || degree>terms+1)
                    throw std::runtime_error("invalid checkpoint polynomial dimensions");
                in.block<K>("polynomial-row",std::span<K>(row_data(r),len*width),len,width);
                lengths[r]=len;degrees[r]=degree;
            }
            processed=terms;
        });
    }
    using Progress=std::function<void(std::size_t,std::size_t)>;
    void process_up_to(std::size_t training,const Progress& progress={}) {
        if(training>moments.size())throw std::out_of_range("recurrence training exceeds available moments");
        for(;processed<training;) {
            validation_active=false;
            const auto t=processed;
            std::fill(discrepancy.begin(),discrepancy.end(),K{});
            // Rows read the same immutable basis and write disjoint discrepancies.
            const bool parallel=threads>1 && t*b*width>=32768;
            (void)parallel;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
            for(std::size_t r=0;r<width;++r) {
                auto* out=discrepancy.data()+r*b;
                for(std::size_t d=0;d<length(r) && d<=t;++d) {
                    const auto* row=coefficient(r,d);
                    if(d==t)for(std::size_t j=0;j<b;++j)out[j]+=row[b+j];
                    for(std::size_t i=0;i<b;++i) {
                        const auto scalar=row[i];
                        if(scalar==K{})continue;
                        const auto* sample=moments[t-d].data()+i*b;
                        for(std::size_t j=0;j<b;++j)out[j]+=scalar*sample[j];
                    }
                }
            }
            std::iota(order.begin(),order.end(),0);
            std::sort(order.begin(),order.end(),[&](auto a,auto c){return degrees[a]!=degrees[c] ? degrees[a]<degrees[c] : a<c;});
            std::size_t pivot_count=0;
            for(auto r:order) {
                auto* out=discrepancy.data()+r*b;
                for(std::size_t pi=0;pi<pivot_count;++pi) {
                    const auto& pivot=pivots[pi];
                    const auto p=pivot.row,c=pivot.column;
                    if(out[c]==K{})continue;
                    const auto factor=out[c]*pivot.inverse;
                    for(std::size_t j=0;j<b;++j)out[j]-=factor*discrepancy[p*b+j];
                    if(length(r)<length(p)) {
                        auto* row=row_data(r);
                        const auto extra=(length(p)-length(r))*width;
                        std::move_backward(row,row+length(r)*width,row+length(p)*width);
                        std::fill_n(row,extra,K{});
                        lengths[r]=length(p);
                    }
                    auto* dst=row_data(r)+(length(r)-length(p))*width;
                    const auto* src=row_data(p);
                    for(std::size_t i=0;i<length(p)*width;++i)dst[i]-=factor*src[i];
                }
                std::size_t c=0;while(c<b && out[c]==K{})++c;
                if(c<b)pivots[pivot_count++]={r,c,out[c].inv()};
            }
            for(std::size_t pi=0;pi<pivot_count;++pi) {
                const auto r=pivots[pi].row;
                std::fill_n(row_data(r)+length(r)*width,width,K{});
                ++lengths[r];
                ++degrees[r];
            }
            ++processed;
            if(progress)progress(processed,training);
        }
    }
    // Cheap structural gate before the full history/holdout check. Merely having
    // b eligible rows is NOT a stopping certificate; generator() still validates.
    bool candidate_ready() const {
        std::size_t eligible=0;
        for(std::size_t r=0;r<width;++r) {
            int p_degree=-1,q_degree=-1;
            for(std::size_t d=0;d<length(r);++d)
                for(std::size_t j=0;j<b;++j) {
                    if(coefficient(r,d)[j]!=K{})p_degree=d;
                    if(coefficient(r,d)[b+j]!=K{})q_degree=d;
                }
            if(p_degree>=0 && q_degree<p_degree && std::size_t(p_degree)<=processed/2)
                if(++eligible==b)return true;
        }
        return false;
    }
    bool generator(packed_generator<K>& result,const Progress& progress={}, bool holdout_first=false) {
        const auto training=processed;
        if(!validation_active)std::fill(validation_prefix.begin(),validation_prefix.end(),0);
        validation_active=true;
        if(result.size()!=b || result.max_degree()<processed/2)throw std::length_error("recurrence output capacity exhausted");
        std::size_t count=0;
        std::iota(order.begin(),order.end(),0);
        std::sort(order.begin(),order.end(),[&](auto a,auto c){return degrees[a]!=degrees[c] ? degrees[a]<degrees[c] : a<c;});
        for(auto r:order) {
            int p_degree=-1,q_degree=-1;
            for(std::size_t d=0;d<length(r);++d)
                for(std::size_t j=0;j<b;++j) {
                    if(coefficient(r,d)[j]!=K{}) p_degree=d;
                    if(coefficient(r,d)[b+j]!=K{}) q_degree=d;
                }
            if(p_degree<0 || q_degree>=p_degree || std::size_t(p_degree)>training/2) continue;

            // Online candidates commonly fit training but fail immediately on a
            // fresh term. Reject those before rescanning the entire fitted prefix.
            // Accepted candidates still undergo the unchanged full-history check.
            if(holdout_first) {
                bool fresh_valid=true;
                for(std::size_t t=training;t<moments.size() && fresh_valid;++t) {
                    for(std::size_t j=0;j<b;++j) {
                        K sum{};
                        for(int d=0;d<=p_degree;++d)
                            for(std::size_t i=0;i<b;++i)sum+=coefficient(r,d)[i]*moments[t-d][i*b+j];
                        if(sum!=K{}) {fresh_valid=false;break;}
                    }
                }
                if(!fresh_valid)continue;
            }
            // Check the entire recurrence, including samples held out of training.
            bool valid=true;
            if(progress)progress(0,moments.size()-p_degree);
            // Validate batches independently, reporting only from the caller.
            for(std::size_t first=std::max(std::size_t(p_degree),validation_prefix[r]);first<moments.size() && valid;first+=256) {
                const auto end=std::min(first+256,moments.size());
                const bool parallel=threads>1 && (end-first)*(p_degree+1)*b*b>=32768;
                (void)parallel;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel) reduction(&&:valid)
                for(std::size_t t=first;t<end;++t) {
                    for(std::size_t j=0;j<b;++j) {
                        K sum{};
                        for(int d=0;d<=p_degree;++d)
                            for(std::size_t i=0;i<b;++i) sum+=coefficient(r,d)[i]*moments[t-d][i*b+j];
                        if(sum!=K{}) { valid=false; break; }
                    }
                }
                if(valid)validation_prefix[r]=end;
                if(progress)progress(end-p_degree,moments.size()-p_degree);
            }
            if(!valid) continue;
            std::copy_n(coefficient(r,0),b,constant.data()+count*b);
            std::copy_n(coefficient(r,p_degree),b,leading.data()+count*b);
            selected[count]=r;selected_degrees[count]=p_degree;
            if(++count==b)break;
        }
        validation_active=false;
        if(count!=b || dense_rank<K>(constant,b)!=b || dense_rank<K>(leading,b)!=b)return false;
        for(std::size_t row=0;row<b;++row) {
            result.set_degree(row,selected_degrees[row]);
            for(std::size_t d=0;d<=selected_degrees[row];++d)
                std::copy_n(coefficient(selected[row],d),b,result.coefficient(row,d).data());
        }
        return true;
    }
};

template <Field K>
std::optional<packed_generator<K>> minimal_generator(
    const packed_moments<K>& moments,std::size_t b,std::size_t training) {
    minimal_generator_state<K> state(moments,b);
    packed_generator<K> result(b,training/2);
    state.process_up_to(training);
    if(!state.generator(result))return std::nullopt;
    return std::optional<packed_generator<K>>(std::move(result));
}
} // namespace VectorSpace::block_wiedemann_detail
