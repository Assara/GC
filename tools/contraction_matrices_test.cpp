#include <cassert>
#include <iostream>
#include "GraphHomology/StackedContraction.hpp"
#include "GraphHomology/NaturalComposition.hpp"
#include "GraphHomology/NaturalRepresentatives.hpp"
#include "VectorSpace/block_wiedemann.hpp"
#include "VectorSpace/sparse_rank.hpp"
using namespace GraphHomology;
std::size_t projected_nonzero_terms = 0;

void storage_checks() {
    using VectorSpace::OwnedArray;
    static_assert(!std::is_copy_constructible_v<OwnedArray<int>>);
    static_assert(std::is_nothrow_move_constructible_v<OwnedArray<int>>);
    OwnedArray<int> empty;
    assert(empty.empty() && empty.begin() == empty.end());
    OwnedArray<int> values(3);
    assert(values[0] == 0);
    values[2] = -7;
    auto moved = std::move(values);
    assert(values.empty() && values.data() == nullptr && moved.back() == -7);
    empty = std::move(moved);
    assert(moved.empty() && empty.size() == 3 && empty.back() == -7);

    ContractionMatrixBuilder m;
    m.rows = 3;
    std::vector<std::pair<std::uint32_t, int>> entries{{2,3},{0,1},{2,-3},{1,-2}};
    m.append_column(entries);
    entries.clear(); m.append_column(entries);
    assert(m.nonzeros() == 2 && m.coefficients[0] == 1 && m.coefficients[1] == -2);
    auto y = m.apply<fieldType>(std::vector<fieldType>{3,7});
    assert(y == std::vector<fieldType>({3,-6,0}));
    entries = {{2,127}}; m.append_column(entries);
    assert(m.coefficients.back() == 127);
    entries = {{1,-127}}; m.append_column(entries);
    assert(m.coefficients.back() == -127);
    auto fixed = std::move(m).finish();
    assert(fixed.nonzeros() == 4 && fixed.columns() == 4);
    assert(fixed.allocated_bytes() == 5 * sizeof(std::size_t) + 4 * 5);
    assert(fixed.apply<fieldType>(std::vector<fieldType>{3,7,0,0}) == y);
    ContractionMatrixBuilder a,b; a.rows = b.rows = 1;
    entries = {{0,1}}; a.append_column(entries); b.append_column(entries);
    bool failed = false;
    try { check_chain(std::move(a).finish(),std::move(b).finish()); } catch (const std::runtime_error&) { failed = true; }
    assert(failed);
}

template <int L, int V, Parity P>
void check_differential(const EnumeratedBasis<L,V,P>& source,
                        const EnumeratedBasis<L,V-1,P>& target, const ContractionMatrix& matrix) {
    for (std::size_t c = 0; c < source.size(); ++c) {
        const auto expected = source.graphs[c].contraction_differential(fieldType{1});
        std::vector<fieldType> values(target.size());
        for (const auto& term : expected) {
            if (term.getCoefficient() == fieldType{0}) continue;
            // Independent oracle: contract and canonicalize first, then inspect
            // the actual child for parallel edges rather than parent triangles.
            const auto& child = term.getValue();
            bool parallel = false;
            for (Int e = 0; e < child.N_EDGES_; ++e)
                for (Int f = 0; f < e; ++f)
                    parallel |= child.getEdge(e) == child.getEdge(f);
            if (parallel) { ++projected_nonzero_terms; continue; }
            values[target.index(child)] += term.getCoefficient();
        }
        for (auto i = matrix.offsets[c]; i < matrix.offsets[c+1]; ++i) {
            assert(values[matrix.row_indices[i]] == fieldType(int(matrix.coefficients[i])));
            values[matrix.row_indices[i]] = 0;
        }
        for (auto x : values) assert(x == fieldType{0});
    }
}

void check_rank(const ContractionMatrix& matrix) {
    compressed_sparse_matrix<fieldType> oracle(matrix.rows);
    for (std::size_t c = 0; c < matrix.columns(); ++c) {
        std::vector<decltype(oracle)::Basis> column;
        for (auto i = matrix.offsets[c]; i < matrix.offsets[c+1]; ++i)
            column.emplace_back(matrix.row_indices[i], fieldType(int(matrix.coefficients[i])));
        oracle.add_col(column);
    }
    const auto expected = VectorSpace::sparse_rank<fieldType>::compute(oracle);
    using Solver = VectorSpace::block_wiedemann_solver<fieldType>;
    const auto actual = Solver(matrix, Solver::options{4,3,8,17}).rank();
    assert(actual.rank == expected.rank && actual.nullity == expected.nullity);
    const auto kernel = Solver(matrix, Solver::options{4,3,8,17}).nullspace();
    assert(kernel.complete && kernel.basis.size() == expected.nullity);
    compressed_sparse_matrix<fieldType> basis(matrix.columns());
    for (const auto& x : kernel.basis) {
        const auto residual = matrix.apply<fieldType>(x);
        assert(std::ranges::all_of(residual, [](auto value) { return value == fieldType{}; }));
        std::vector<decltype(basis)::Basis> column;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (x[i] != fieldType{}) column.emplace_back(i, x[i]);
        basis.add_col(column);
    }
    assert(VectorSpace::sparse_rank<fieldType>::compute(basis).rank == expected.nullity);
}

void check_stack(const ContractionMatrix& down, const ContractionMatrix& up) {
    using K = fieldType;
    using Solver = VectorSpace::block_wiedemann_solver<K>;
    StackedContraction stack(down, up);
    // Independently assemble a small dense oracle from the two sparse inputs.
    std::vector<std::vector<K>> dense(stack.rows(), std::vector<K>(stack.columns()));
    for (std::size_t c=0;c<down.columns();++c)
        for (auto i=down.offsets[c];i<down.offsets[c+1];++i)
            dense[down.row_indices[i]][c] += K(int(down.coefficients[i]));
    for (std::size_t c=0;c<up.columns();++c)
        for (auto i=up.offsets[c];i<up.offsets[c+1];++i)
            dense[down.rows+c][up.row_indices[i]] += K(int(up.coefficients[i]));
    constexpr std::size_t b=3;
    std::vector<K> x(stack.columns()*b), y(stack.rows()*b), ax(y.size()), aty(x.size());
    for (std::size_t i=0;i<x.size();++i) x[i]=K(i+1);
    for (std::size_t i=0;i<y.size();++i) y[i]=K(i+3);
    stack.apply<K>(x,ax,b); stack.transpose<K>(y,aty,b);
    for (std::size_t r=0;r<stack.rows();++r) for (std::size_t j=0;j<b;++j) {
        K expected{};
        for (std::size_t c=0;c<stack.columns();++c) expected+=dense[r][c]*x[c*b+j];
        assert(ax[r*b+j]==expected);
    }
    compressed_sparse_matrix<K> oracle(stack.rows());
    for (std::size_t c=0;c<stack.columns();++c) {
        for (std::size_t j=0;j<b;++j) {
            K expected{};
            for (std::size_t r=0;r<stack.rows();++r) expected+=dense[r][c]*y[r*b+j];
            assert(aty[c*b+j]==expected);
        }
        std::vector<decltype(oracle)::Basis> column;
        for (std::size_t r=0;r<stack.rows();++r)
            if (dense[r][c]!=K{}) column.emplace_back(r,dense[r][c]);
        oracle.add_col(column);
    }
    auto kernel=Solver(stack.rows(),stack.columns(),
        [&stack](auto in,auto out,auto width){stack.apply<K>(in,out,width);},
        [&stack](auto in,auto out,auto width){stack.transpose<K>(in,out,width);}).nullspace();
    assert(kernel.complete && kernel.basis.size()==VectorSpace::sparse_rank<K>::compute(oracle).nullity);
    compressed_sparse_matrix<K> basis(stack.columns());
    for (const auto& z:kernel.basis) {
        for (const auto& row:dense) {
            K residual{};
            for (std::size_t c=0;c<z.size();++c) residual+=row[c]*z[c];
            assert(residual==K{});
        }
        std::vector<decltype(basis)::Basis> column;
        for (std::size_t c=0;c<z.size();++c) if(z[c]!=K{})column.emplace_back(c,z[c]);
        basis.add_col(column);
    }
    assert(VectorSpace::sparse_rank<K>::compute(basis).rank==kernel.basis.size());
}

// Independent sign/multiplicity oracle: enumerate raw splits only in tests.
template<int L, int V, Parity P>
void check_natural_adjoint(const EnumeratedBasis<L,V,P>& source,
                           const EnumeratedBasis<L,V+1,P>& target,
                           const NaturalAdjointMatrix& matrix) {
    using G = typename EnumeratedBasis<L,V,P>::G;
    for (std::size_t c=0;c<source.size();++c) {
        std::vector<fieldType> expected(target.size());
        const auto splits = source.graphs[c].split_vertex_differential(fieldType{1});
        for (const auto& term : splits) {
            if (term.getCoefficient() == fieldType{}) continue;
            if (GraphGeneration::CutVertexSplitRule<V+1>(term.getValue()).vertex >= 0) continue;
            expected[target.index(term.getValue())] += term.getCoefficient();
        }
        for(auto i=matrix.offsets[c];i<matrix.offsets[c+1];++i) {
            if (expected[matrix.row_indices[i]] != fieldType(matrix.coefficients[i])) {
                std::cerr << "adjoint mismatch parity=" << int(P) << " V=" << V
                    << " column=" << c << " row=" << matrix.row_indices[i]
                    << " expected=" << expected[matrix.row_indices[i]].value()
                    << " actual=" << matrix.coefficients[i] << std::endl;
                std::abort();
            }
            expected[matrix.row_indices[i]] = fieldType{};
        }
        for(auto x:expected) assert(x==fieldType{});
    }
}

void check_natural_composition(const ContractionMatrix& down,const ContractionMatrix& up,
                               const NaturalAdjoints& adjoints) {
    for(std::size_t b:{1,8}) {
        NaturalComposition<fieldType> op32(down,up,adjoints,b,1);
        NaturalComposition<fieldType,std::int64_t> op64(down,up,adjoints,b,8);
        std::vector<fieldType> x(down.columns()*b), actual(x.size()),wide(x.size()),expected(x.size());
        for(std::size_t i=0;i<x.size();++i) x[i]=fieldType(i%2 ? -1 : i*739);
        op32.apply(x,actual); op64.apply(x,wide);
        // Field arithmetic after every addition, with independent column loops.
        const auto transposed = [b](const auto& m,const auto& in) {
            std::vector<fieldType> out(m.columns()*b);
            for(std::size_t c=0;c<m.columns();++c)
                for(auto i=m.offsets[c];i<m.offsets[c+1];++i)
                    for(std::size_t j=0;j<b;++j)
                        out[c*b+j] += fieldType(m.coefficients[i])*in[m.row_indices[i]*b+j];
            return out;
        };
        const auto low=transposed(down,transposed(adjoints.down,x));
        const auto high=transposed(adjoints.up,transposed(up,x));
        for(std::size_t i=0;i<x.size();++i) expected[i]=low[i]+high[i];
        assert(actual==expected && wide==expected);
    }
}

void check_natural_solver(const ContractionMatrix& down,const ContractionMatrix& up,const NaturalAdjoints& adjoints) {
    using K=fieldType;
    using Solver=VectorSpace::block_wiedemann_solver<K>;
    NaturalComposition<K> composition(down,up,adjoints,8);
    auto solver=Solver::from_square_operator(down.columns(),[&](auto in,auto out,auto b){
        composition.apply(in,out,b);
    });
    auto kernel=solver.nullspace();
    assert(kernel.complete);
    convert_natural_representatives<K>(kernel.basis,down,adjoints);
    compressed_sparse_matrix<K> image(down.columns());
    for(std::size_t c=0;c<up.columns();++c) {
        std::vector<decltype(image)::Basis> col;
        for(auto i=up.offsets[c];i<up.offsets[c+1];++i)col.emplace_back(up.row_indices[i],K(up.coefficients[i]));
        image.add_col(col);
    }
    const auto image_rank=VectorSpace::sparse_rank<K>::compute(image).rank;
    for(const auto& x:kernel.basis) {
        std::vector<decltype(image)::Basis> col;
        for(std::size_t i=0;i<x.size();++i)if(x[i]!=K{})col.emplace_back(i,x[i]);
        image.add_col(col);
    }
    assert(VectorSpace::sparse_rank<K>::compute(image).rank==image_rank+kernel.basis.size());
    compressed_sparse_matrix<K> differential(down.rows);
    for(std::size_t c=0;c<down.columns();++c) {
        std::vector<decltype(differential)::Basis> col;
        for(auto i=down.offsets[c];i<down.offsets[c+1];++i)col.emplace_back(down.row_indices[i],K(down.coefficients[i]));
        differential.add_col(col);
    }
    assert(kernel.basis.size()==VectorSpace::sparse_rank<K>::compute(differential).nullity-image_rank);
}

void conversion_checks() {
    ContractionMatrixBuilder d,u;
    d.rows=1;u.rows=3;
    std::vector<std::pair<std::uint32_t,int>> entry{{0,1}};
    for(int i=0;i<3;++i)d.append_column(entry);
    entry={{0,1},{1,-1}};u.append_column(entry);
    auto down=std::move(d).finish(),up=std::move(u).finish();
    NaturalAdjoints adj;
    adj.lower_aut=AutomorphismSizes(std::span<const std::uint64_t>(std::array<std::uint64_t,1>{6}));
    adj.middle_aut=AutomorphismSizes(std::span<const std::uint64_t>(std::array<std::uint64_t,3>{1,2,3}));
    adj.upper_aut=AutomorphismSizes(std::span<const std::uint64_t>(std::array<std::uint64_t,1>{1}));
    adj.down=natural_adjoint(down,adj.middle_aut,adj.lower_aut);
    adj.up=natural_adjoint(up,adj.upper_aut,adj.middle_aut);
    check_natural_solver(down,up,adj);
    std::vector<std::vector<fieldType>> basis{{2,2,-9}};
    convert_natural_representatives<fieldType>(basis,down,adj);
    assert(basis[0]==std::vector<fieldType>({1,fieldType{2}.inv(),fieldType{-3}/fieldType{2}}));
    bool rejected=false;
    basis={{1,0,0}};
    try {convert_natural_representatives<fieldType>(basis,down,adj);}
    catch(const std::runtime_error&) {rejected=true;}
    assert(rejected);
    adj.middle_aut[0]=fieldType::characteristic();
    rejected=false;
    try {check_automorphism_units<fieldType>(adj);}
    catch(const std::runtime_error&) {rejected=true;}
    assert(rejected);
}

void natural_safety_checks() {
    ContractionMatrixBuilder builder;
    builder.rows=1;
    std::vector<std::pair<std::uint32_t,int>> entries{{0,-1}};
    builder.append_column(entries);
    auto down=std::move(builder).finish();
    ContractionMatrix up;
    up.rows=1;
    NaturalAdjoints adj;
    const std::array<std::uint64_t,1> source{1},target{70000},two{2};
    adj.down=natural_adjoint(down,source,target);
    adj.up=natural_adjoint(up,std::span<const std::uint64_t>{},source);
    bool rejected=false;
    try { NaturalComposition<fieldType,std::int32_t> narrow(down,up,adj,1); }
    catch(const std::overflow_error&) { rejected=true; }
    assert(rejected);
    NaturalComposition<fieldType,std::int64_t> wide(down,up,adj,1);
    std::vector<fieldType> input{fieldType{-1}},output(1);
    wide.apply(input,output);
    assert(output[0]==fieldType{-70000});
    rejected=false;
    try { auto invalid=natural_adjoint(down,two,source); }
    catch(const std::runtime_error&) { rejected=true; }
    assert(rejected);
    const std::array<std::uint64_t,1> large{UINT64_MAX};
    auto cancelled=natural_adjoint(down,large,large);
    assert(cancelled.coefficients[0]==-1);
    rejected=false;
    try { auto invalid=natural_adjoint(down,source,large); }
    catch(const std::overflow_error&) { rejected=true; }
    assert(rejected);
}

template <Parity P, int V = 5>
void windows(const char* directory) {
    ContractionWindow<6,V,P> w(directory);
    check_differential(w.middle,w.lower,w.down);
    check_differential(w.upper,w.middle,w.up);
    check_chain(w.down,w.up);
    check_rank(w.down); check_rank(w.up);
    check_stack(w.down,w.up);
    NaturalAdjoints adjoints(w);
    if constexpr (V == 8) {
        // Brute-force vertex permutations independently check automorphism counts.
        for(std::size_t i=0;i<w.middle.size();++i) {
            bool adjacent[V][V]{};
            for(int e=0;e<decltype(w.middle)::G::N_EDGES_;++e) {
                auto [a,b]=w.middle.graphs[i].getEdge(e);
                adjacent[a][b]=adjacent[b][a]=true;
            }
            std::array<int,V> permutation;
            std::iota(permutation.begin(),permutation.end(),0);
            std::uint64_t count=0;
            do {
                bool same=true;
                for(int a=0;a<V && same;++a)
                    for(int b=0;b<a;++b)
                        if(adjacent[a][b]!=adjacent[permutation[a]][permutation[b]]) { same=false; break; }
                count+=same;
            } while(std::next_permutation(permutation.begin(),permutation.end()));
            assert(count==adjoints.middle_aut[i]);
        }
    }
    check_natural_adjoint(w.lower,w.middle,adjoints.down);
    check_natural_adjoint(w.middle,w.upper,adjoints.up);
    check_natural_composition(w.down,w.up,adjoints);
    check_natural_solver(w.down,w.up,adjoints);
    std::cout << "L6 " << (P == Parity::even ? "even" : "odd") << " V=" << V << ": matrices, ranks, nullspace bases and unweighted stack verified; d^2=0\n";
    if constexpr (V < 10) windows<P,V+1>(directory);
}
void parallel_stack_checks() {
    ContractionMatrixBuilder down_builder,up_builder;
    down_builder.rows=513;up_builder.rows=1024;
    auto fill=[](ContractionMatrixBuilder& m,std::size_t columns) {
        for(std::size_t c=0;c<columns;++c) {
            std::vector<std::pair<std::uint32_t,int>> entries;
            for(std::size_t j=0;j<11;++j)
                entries.emplace_back((c*17+j*31)%m.rows,int((c+j)%7)-3);
            m.append_column(entries);
        }
    };
    fill(down_builder,1024);fill(up_builder,731);
    auto down = std::move(down_builder).finish();
    auto up = std::move(up_builder).finish();
    NaturalAdjoints natural;
    std::vector<std::uint64_t> lower(down.rows,1),middle(down.columns(),1),upper(up.columns(),1);
    natural.down=natural_adjoint(down,middle,lower);
    natural.up=natural_adjoint(up,upper,middle);
    check_natural_composition(down,up,natural);
    StackedContraction serial(down,up,1),parallel(down,up,8);
    for(std::size_t b:{1,8}) {
        std::vector<fieldType> x(serial.columns()*b),y(serial.rows()*b);
        for(std::size_t i=0;i<x.size();++i)x[i]=fieldType(i+1);
        for(std::size_t i=0;i<y.size();++i)y[i]=fieldType(3*i+7);
        std::vector<fieldType> ax(y.size()),px(y.size()),ty(x.size()),pty(x.size());
        serial.apply<fieldType>(x,ax,b);parallel.apply<fieldType>(x,px,b);
        serial.transpose<fieldType>(y,ty,b);parallel.transpose<fieldType>(y,pty,b);
        assert(ax==px && ty==pty);
        // Independent scatter oracle, including both contributions per output.
        std::vector<fieldType> expected_ax(y.size()),expected_ty(x.size());
        for(std::size_t c=0;c<down.columns();++c)
            for(auto i=down.offsets[c];i<down.offsets[c+1];++i)
                for(std::size_t j=0;j<b;++j) {
                    expected_ax[down.row_indices[i]*b+j]+=x[c*b+j]*down.coefficients[i];
                    expected_ty[c*b+j]+=y[down.row_indices[i]*b+j]*down.coefficients[i];
                }
        for(std::size_t c=0;c<up.columns();++c)
            for(auto i=up.offsets[c];i<up.offsets[c+1];++i)
                for(std::size_t j=0;j<b;++j) {
                    expected_ax[(down.rows+c)*b+j]+=x[up.row_indices[i]*b+j]*up.coefficients[i];
                    expected_ty[up.row_indices[i]*b+j]+=y[(down.rows+c)*b+j]*up.coefficients[i];
                }
        assert(ax==expected_ax && ty==expected_ty);
    }
    std::cout << "Stored transpose: 1/8-thread products match independent scatter oracle\n";
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    storage_checks();
    natural_safety_checks();
    conversion_checks();
    parallel_stack_checks();
    windows<Parity::even>(argv[1]);
    windows<Parity::odd>(argv[1]);
    assert(projected_nonzero_terms > 0);
    std::cout << "Natural adjoints match direct splitting; automorphism, int32/int64 composition and overflow checks passed\n";
    std::cout << "Projected out " << projected_nonzero_terms << " nonzero parallel-edge terms\n";
}
