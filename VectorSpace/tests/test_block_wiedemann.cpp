#include <iostream>
#include <random>
#include <stdexcept>
#include "VectorSpace/block_wiedemann.hpp"
#include "VectorSpace/sparse_rank.hpp"
using K=fieldType;
using Matrix=compressed_sparse_matrix<K>;
using Solver=VectorSpace::block_wiedemann_solver<K>;
std::size_t cases=0;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
Matrix matrix(const std::vector<std::vector<K>>& a,std::size_t n) {
    Matrix m(a.size());
    for(std::size_t c=0;c<n;++c) {
        std::vector<Matrix::Basis> col;
        for(std::size_t r=0;r<a.size();++r)if(a[r][c]!=K{})col.emplace_back(r,a[r][c]);
        m.add_col(col);
    }
    return m;
}
void check(const Matrix& m,std::size_t block,std::uint64_t seed=17) {
    Solver::options config{block,3,8,seed};
    config.sequence_capacity=2*std::max(m.image_dim(),m.domain_dim())+10;
    Solver solver(m,config);
    auto rank=solver.rank();
    auto exact=VectorSpace::sparse_rank<K>::compute(m);
    if(rank.rank!=exact.rank){std::cerr<<"rows="<<m.image_dim()<<" cols="<<m.domain_dim()<<" block="<<block<<" got="<<rank.rank<<" expected="<<exact.rank<<'\n';throw std::runtime_error("rank mismatch");}
    require(rank.nullity==exact.nullity,"nullity mismatch");
    auto kernel=solver.nullspace();
    require(kernel.complete,"nullspace extraction incomplete");
    require(kernel.basis.size()==exact.nullity,"nullspace dimension mismatch");
    Matrix basis(m.domain_dim());
    for(const auto& x:kernel.basis) {
        require(x.size()==m.domain_dim(),"nullspace coordinate count");
        std::vector<K> actual(m.image_dim());
        std::vector<Matrix::Basis> column;
        for(std::size_t c=0;c<m.domain_dim();++c) {
            if(x[c]!=K{})column.emplace_back(c,x[c]);
            for(const auto& t:m.get_column(c))actual[t.getValue()]+=x[c]*t.getCoefficient();
        }
        for(auto value:actual)require(value==K{},"invalid nullspace vector");
        basis.add_col(column);
    }
    require(VectorSpace::sparse_rank<K>::compute(basis).rank==kernel.basis.size(),"dependent nullspace vectors");
    auto rhs=std::make_unique<K[]>(m.image_dim());
    for(std::size_t c=0;c<m.domain_dim();++c)
        for(auto& x:m.get_column(c))rhs[x.getValue()]+=K(c+1)*x.getCoefficient();
    auto solutions=solver.solve_MX_equals_y(rhs);
    require(!solutions.empty(),"consistent system was not solved");
    for(auto& x:solutions) {
        std::vector<K> actual(m.image_dim());
        for(std::size_t c=0;c<m.domain_dim();++c)
            for(auto& t:m.get_column(c))actual[t.getValue()]+=x[c]*t.getCoefficient();
        for(std::size_t r=0;r<actual.size();++r)require(actual[r]==rhs[r],"invalid returned solution");
    }
    ++cases;
}
void check_direct_square() {
    const std::vector<std::vector<std::vector<K>>> examples{
        {{1,1},{2,2}}, {{0,1},{1,0}}, {{0,0},{0,0}},
        {{1,2,3},{2,4,6},{3,6,9}}
    };
    for(const auto& dense:examples) for(std::size_t b:{1,2,8}) for(bool online:{false,true}) {
        const auto n=dense.size();
        const auto exact=VectorSpace::sparse_rank<K>::compute(matrix(dense,n));
        Solver::options config{b,1,8,17};
        config.incremental_recurrence=online;config.recurrence_check_interval=2;
        auto solver=Solver::from_square_operator(n,[&](auto in,auto out,std::size_t width){
            std::fill(out.begin(),out.end(),K{});
            for(std::size_t r=0;r<n;++r)for(std::size_t c=0;c<n;++c)
                for(std::size_t j=0;j<width;++j)out[r*width+j]+=dense[r][c]*in[c*width+j];
        },config);
        const auto kernel=solver.nullspace();
        require(kernel.complete && kernel.rank_estimate.rank==exact.rank,"direct square rank/nullity");
        Matrix independent(n);
        for(const auto& x:kernel.basis) {
            for(const auto& row:dense) {
                K residual{};
                for(std::size_t c=0;c<n;++c)residual+=row[c]*x[c];
                require(residual==K{},"direct square residual");
            }
            std::vector<Matrix::Basis> column;
            for(std::size_t c=0;c<n;++c)if(x[c]!=K{})column.emplace_back(c,x[c]);
            independent.add_col(column);
        }
        require(VectorSpace::sparse_rank<K>::compute(independent).rank==exact.nullity,"direct square independence");
        ++cases;
    }
}
void check_early_stop() {
    std::size_t steps[2]{};
    for(int mode=0;mode<2;++mode) {
        Solver::options config;
        config.incremental_recurrence=mode;
        auto solver=Solver::from_square_operator(128,[](auto in,auto out,std::size_t b){
            for(std::size_t i=0;i<128;++i)for(std::size_t j=0;j<b;++j)
                out[i*b+j]=i<8 ? K(i+1)*in[i*b+j] : K{};
        },config);
        const auto result=solver.rank();
        require(result.rank==8,"early stopping exact diagonal rank");
        require(solver.sequence_stats.moments==solver.sequence_stats.recurrence_updates+config.holdout,
            "holdout terms entered training");
        steps[mode]=solver.sequence_stats.moments;
    }
    require(steps[1]<steps[0],"incremental did not stop early");
    ++cases;
}
int main() {
    try {
        check_direct_square();
        check_early_stop();
        std::mt19937_64 rng(41);
        // Both moment layouts, uneven widths, and the real L9 domain size.
        for(std::size_t n:{0,13,68987})for(std::size_t b:{1,3,8}) {
            std::vector<K> u(n*b),v(n*b),expected(b*b);
            for(auto& x:u)x=K::sample(rng);
            for(auto& x:v)x=K::sample(rng);
            for(std::size_t r=0;r<n;++r)
                for(std::size_t i=0;i<b;++i)for(std::size_t j=0;j<b;++j)
                    expected[i*b+j]+=u[r*b+i]*v[r*b+j];
            for(bool transposed:{false,true})for(int threads:{1,8}) {
                std::vector<K> actual(b*b,K(7));
                VectorSpace::block_wiedemann_detail::project_block<K>(u,v,actual,n,b,transposed,threads);
                for(std::size_t i=0;i<b;++i)for(std::size_t j=0;j<b;++j)
                    require(actual[transposed?j*b+i:i*b+j]==expected[i*b+j],"block projection mismatch");
                if constexpr(std::same_as<K,Z32783>) {
                    std::vector<K> pu(n*b),pv(n*b);
                    VectorSpace::block_wiedemann_detail::pack_projection<K>(u,pu,n,b,threads);
                    VectorSpace::block_wiedemann_detail::pack_projection<K>(v,pv,n,b,threads);
                    for(std::size_t lane=0;lane<b;++lane)for(std::size_t r=0;r<n;++r) {
                        require(pu[lane*n+r]==u[r*b+lane],"projection packing mismatch");
                        require(pv[lane*n+r]==v[r*b+lane],"current packing mismatch");
                    }
                    VectorSpace::block_wiedemann_detail::project_packed32783(pu,pv,actual,n,b,transposed,threads);
                    for(std::size_t i=0;i<b;++i)for(std::size_t j=0;j<b;++j)
                        require(actual[transposed?j*b+i:i*b+j]==expected[i*b+j],"packed projection mismatch");
                }
            }
        }

        for(std::size_t rows:{0,1,3,8,13,24})for(std::size_t cols:{0,1,5,12,19}) {
            for(std::size_t rank:{std::size_t{0},std::min(rows,cols)/2,std::min(rows,cols)}) {
                std::vector<std::vector<K>> l(rows,std::vector<K>(rank)),r(rank,std::vector<K>(cols));
                for(auto& row:l)for(auto& x:row)x=K::sample(rng);
                for(auto& row:r)for(auto& x:row)x=K::sample(rng);
                std::vector<std::vector<K>> a(rows,std::vector<K>(cols));
                for(std::size_t i=0;i<rows;++i)for(std::size_t j=0;j<cols;++j)
                    for(std::size_t k=0;k<rank;++k)a[i][j]+=l[i][k]*r[k][j];
                auto m=matrix(a,cols);
                for(auto b:{1,2,4})check(m,b);
            }
        }
        // A nonzero row with A A^T = 0: rank cannot be obtained from plain Gram.
        std::vector<int> roots(K::characteristic(),-1);
        for(std::size_t i=0;i<roots.size();++i)roots[(K(i)*K(i)).value()]=i;
        bool found=false;
        for(std::size_t i=0;i<roots.size() && !found;++i) {
            auto j=roots[(-K{1}-K(i)*K(i)).value()];
            if(j<0)continue;
            auto m=matrix({{K{1},K(i),K(j)}},3);
            check(m,2);check(matrix({{K{1}},{K(i)},{K(j)}},1),2);
            found=true;
        }
        require(found,"isotropic test setup failed");
        // Legacy LIL entry point remains usable.
        lil_matrix<K> legacy;legacy.add_element(0,0,1);legacy.add_element(1,1,1);
        Solver legacy_solver(legacy,2);
        require(legacy_solver.rank().rank==2,"legacy constructor rank");
        auto no=matrix({{K{1}},{K{0}}},1);
        auto rhs=std::make_unique<K[]>(2);rhs[1]=1;
        require(Solver(no).solve_MX_equals_y(rhs).empty(),"inconsistent system returned a solution");
        bool rejected=false;
        try{Solver invalid(no,Solver::options{0,3,8,17});}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"zero block size accepted");
        // Square non-symmetric operator: cached reconstruction must use A^T A,
        // not A A^T, even though the two dimensions are equal.
        auto square=matrix({{1,2,3},{0,1,1},{1,3,4}},3);
        check(square,2);
        Solver square_solver(square,Solver::options{2,3,8,17});
        auto square_kernel=square_solver.nullspace();
        require(square_kernel.complete && square_kernel.basis.size()==1,"square cached kernel");
        require(square_solver.sequence_stats.sequences==3,"square cache was not reused");
        // Incremental recurrence state must survive extension of the moment array.
        VectorSpace::block_wiedemann_detail::packed_moments<K> moments(2,24);
        VectorSpace::block_wiedemann_detail::minimal_generator_state<K> state(moments,2);
        for(std::size_t t=0;t<24;++t) {
            std::vector<K> moment(4);
            for(std::size_t a=0;a<7;++a) {
                K power=1;for(std::size_t i=0;i<t;++i)power*=K(a+1);
                for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j)
                    moment[i*2+j]+=K((a+1)*(i+1)+3)*power*K((a+2)*(j+1)+1);
            }
            std::ranges::copy(moment,moments.append().begin());
            if(t==5)state.process_up_to(4);
            if(t==12)state.process_up_to(10);
        }
        state.process_up_to(16);
        VectorSpace::block_wiedemann_detail::packed_generator<K> incremental(2,12);
        const bool incremental_ok=state.generator(incremental);
        auto batch=VectorSpace::block_wiedemann_detail::minimal_generator(moments,2,16);
        require(incremental_ok && batch.has_value(),"incremental recurrence missing");
        for(std::size_t i=0;i<2;++i)
            require(std::ranges::equal(incremental.row(i),batch->row(i)),"incremental recurrence changed");
        // Exercise parallel recurrence construction and held-out validation.
        {
            constexpr std::size_t b=8,degree=128,training=2*degree+2;
            VectorSpace::block_wiedemann_detail::packed_moments<K> samples(b,training+8);
            for(std::size_t t=0;t<training+8;++t)samples.append();
            for(std::size_t lane=0;lane<b;++lane)
                for(std::size_t a=0;a<degree;++a) {
                    K power=1,lambda=K(1+lane*degree+a);
                    for(std::size_t t=0;t<samples.size();++t) {samples[t][lane*b+lane]+=power;power*=lambda;}
                }
            using State=VectorSpace::block_wiedemann_detail::minimal_generator_state<K>;
            State serial(samples,b,1),parallel(samples,b,8);
            serial.process_up_to(training);
            parallel.process_up_to(training/2);
            parallel.process_up_to(training);
            VectorSpace::block_wiedemann_detail::packed_generator<K> expected(b,training/2),actual(b,training/2),online(b,training/2);
            const bool expected_ok=serial.generator(expected),actual_ok=parallel.generator(actual);
            require(expected_ok && actual_ok,"parallel recurrence missing");
            for(std::size_t i=0;i<b;++i)
                require(std::ranges::equal(expected.row(i),actual.row(i)),"parallel recurrence differs");
            const bool online_ok=parallel.generator(online,{},true);
            require(online_ok,"holdout-first valid recurrence missing");
            for(std::size_t i=0;i<b;++i)
                require(std::ranges::equal(expected.row(i),online.row(i)),"holdout-first recurrence differs");
            samples[samples.size()-1][0]+=K(1);
            require(!serial.generator(expected) && !parallel.generator(actual),"corrupt holdout accepted");
            require(!serial.generator(expected,{},true) && !parallel.generator(actual,{},true),"holdout-first corrupt sequence accepted");
        }
        // Every signed-byte coefficient, including values beyond +/-1.
        for(int scalar=-128;scalar<=127;++scalar) {
            K in[3]{0,17,-5},actual[3]{3,4,5},expected[3]{3,4,5};
            for(std::size_t i=0;i<3;++i)expected[i]+=in[i]*K(scalar);
            VectorSpace::add_signed_block(actual,in,3,SmallSignedInt(scalar));
            for(std::size_t i=0;i<3;++i)require(actual[i]==expected[i],"signed block multiply");
        }
        // Enough work to exercise the OpenMP field-valued gather path.
        check(matrix(std::vector<std::vector<K>>(129,std::vector<K>(65,K{1})),65),8);
        // Native int8 CSC and a matrix-free vertical stack share the API.
        struct ByteMatrix {
            std::size_t rows=2;
            std::vector<std::size_t> offsets{0,1,2,4};
            std::vector<std::uint32_t> row_indices{0,1,0,1};
            std::vector<SmallSignedInt> coefficients{1,-1,1,1};
            std::size_t columns() const{return offsets.size()-1;}
        } bytes;
        auto native=Solver(bytes).nullspace();
        require(native.complete && native.basis.size()==1,"native byte nullspace");
        const auto& x=native.basis.front();
        require(x[0]+x[2]==K{} && -x[1]+x[2]==K{},"native byte residual");
        ByteMatrix large_bytes;
        large_bytes.rows=129;
        large_bytes.offsets.assign(1,0);large_bytes.row_indices.clear();large_bytes.coefficients.clear();
        for(std::size_t c=0;c<65;++c) {
            for(std::size_t r=0;r<129;++r) {
                large_bytes.row_indices.push_back(r);large_bytes.coefficients.push_back(1);
            }
            large_bytes.offsets.push_back(large_bytes.coefficients.size());
        }
        auto large_kernel=Solver(large_bytes).nullspace();
        require(large_kernel.complete && large_kernel.basis.size()==64,"parallel native kernel");
        for(const auto& z:large_kernel.basis) {
            K sum{};for(auto value:z)sum+=value;
            require(sum==K{},"parallel native residual");
        }
        // [A; B^T] = [[1,1,0]; [0,0,1]]; no weighting of either block.
        Solver stacked(2,3,[](auto in,auto out,std::size_t b){
            for(std::size_t j=0;j<b;++j){out[j]=in[j]+in[b+j];out[b+j]=in[2*b+j];}
        },[](auto in,auto out,std::size_t b){
            for(std::size_t j=0;j<b;++j){out[j]=out[b+j]=in[j];out[2*b+j]=in[b+j];}
        });
        auto stack_kernel=stacked.nullspace();
        require(stack_kernel.complete && stack_kernel.basis.size()==1,"stacked nullspace");
        require(stack_kernel.basis[0][0]+stack_kernel.basis[0][1]==K{} && stack_kernel.basis[0][2]==K{},"stacked residual");
        std::cout<<"Block Wiedemann: "<<cases<<" rank, nullspace and verified-solution cases passed; byte CSC and stack passed\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
