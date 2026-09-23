#include <chrono>
#include <cstdlib>
#include <iostream>
#include "GraphHomology/StackedContraction.hpp"
#include "GraphHomology/NaturalComposition.hpp"
#include "GraphHomology/NaturalRepresentatives.hpp"
#include "VectorSpace/block_wiedemann.hpp"
#ifndef GC_HOMOLOGY_LOOP
#define GC_HOMOLOGY_LOOP 6
#endif
#ifndef GC_HOMOLOGY_VERTICES
#define GC_HOMOLOGY_VERTICES 8
#endif

template <GraphHomology::Parity P>
void run(const char* directory, bool ranks, std::size_t block_size, const char* output, std::uint64_t seed, bool natural, bool legacy, bool prepare) {
    constexpr int L = GC_HOMOLOGY_LOOP, V = GC_HOMOLOGY_VERTICES;
    const auto start = std::chrono::steady_clock::now();
    GraphHomology::ContractionWindow<L, V, P> window(directory);
    std::cout << "L=" << L << " V=" << V << " parity=" << (P == GraphHomology::Parity::even ? "even" : "odd")
        << " field=" << fieldType::name() << '\n';
    const auto basis = [](int v, const auto& b) {
        std::cout << "basis V=" << v << " input=" << b.input_graphs << " sign_zero=" << b.zero_graphs
            << " dimension=" << b.size() << '\n';
    };
    basis(V - 1, window.lower); basis(V, window.middle); basis(V + 1, window.upper);
    const auto matrix = [](const char* name, const auto& m) {
        std::cout << name << " rows=" << m.rows << " columns=" << m.columns() << " nnz=" << m.nonzeros()
            << " allocated_bytes=" << m.allocated_bytes() << '\n';
    };
    matrix("down", window.down); matrix("up", window.up);
    std::cout << "total_allocated_bytes=" << window.allocated_bytes() << " build_seconds="
        << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << std::endl;
    GraphHomology::check_chain(window.down, window.up);
    std::cout << "Verified down * up = 0 over the integers\n";
    if (natural) {
        const auto adjoint_start = std::chrono::steady_clock::now();
        GraphHomology::NaturalAdjoints adjoints(window);
        matrix("natural_down_adjoint",adjoints.down);
        matrix("natural_up_adjoint",adjoints.up);
        GraphHomology::NaturalComposition<fieldType> composition(window.down,window.up,adjoints,block_size);
        std::cout << "natural_adjoint_bytes=" << adjoints.allocated_bytes()
            << " build_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-adjoint_start).count()
            << " accumulator_bits=" << 8*sizeof(GraphAccumulator)
            << " path_bound=" << composition.path_bound()
            << " workspace_bytes=" << composition.workspace_bytes() << std::endl;
        std::vector<fieldType> x(window.middle.size()*block_size), y(x.size());
        for(std::size_t i=0;i<x.size();++i) x[i]=fieldType(i+seed);
        const auto multiply_start = std::chrono::steady_clock::now();
        composition.apply(x,y);
        std::cout << "natural_composition_seconds="
            << std::chrono::duration<double>(std::chrono::steady_clock::now()-multiply_start).count() << std::endl;
    }
    if (output) {
        using Solver = VectorSpace::block_wiedemann_solver<fieldType>;
        const auto started = std::chrono::steady_clock::now();
        typename Solver::options solver_options{block_size,1,8,seed,8,&std::cout};
        solver_options.checkpoint_path=std::string(output)+".recurrence";
        if(const auto* sequence=std::getenv("GC_RECURRENCE_SEQUENCE"))solver_options.recurrence_sequence_path=sequence;
        if(const auto* seconds=std::getenv("GC_CHECKPOINT_SECONDS")) {
            std::size_t end=0;const auto interval=std::stod(seconds,&end);
            if(end!=std::string(seconds).size() || !std::isfinite(interval) || interval<0)
                throw std::invalid_argument("GC_CHECKPOINT_SECONDS must be finite and nonnegative");
            solver_options.checkpoint_seconds=interval;
            solver_options.reconstruction_checkpoint_seconds=interval;
        }
        if(const auto* capacity=std::getenv("GC_SEQUENCE_CAPACITY")) {
            const std::string value(capacity);
            if(value.empty() || value.find_first_not_of("0123456789")!=std::string::npos)
                throw std::invalid_argument("GC_SEQUENCE_CAPACITY must be a positive integer");
            solver_options.sequence_capacity=std::stoull(value);
            if(!solver_options.sequence_capacity)throw std::invalid_argument("GC_SEQUENCE_CAPACITY must be positive");
        }
        if(const auto* resume=std::getenv("GC_RESUME_CHECKPOINT")) {
            if(std::string(resume)!="1")throw std::invalid_argument("GC_RESUME_CHECKPOINT must be 1");
            solver_options.resume_checkpoint=true;
        }
        typename Solver::nullspace_result kernel;
        GraphHomology::AutomorphismSizes saved_aut;
        std::size_t operator_rows=window.middle.size();
        if(!legacy) {
            GraphHomology::NaturalAdjoints adjoints(window);
            GraphHomology::check_automorphism_units<fieldType>(adjoints);
            GraphHomology::NaturalComposition<fieldType> composition(window.down,window.up,adjoints,block_size);
            std::cout << "operator=transposed_natural_composition accumulator_bits=" << 8*sizeof(GraphAccumulator)
                << " path_bound=" << composition.path_bound() << " workspace_bytes=" << composition.workspace_bytes()
                << " adjoint_bytes=" << adjoints.allocated_bytes() << " threads=8" << std::endl;
            auto solver=Solver::from_square_operator(window.middle.size(),
                [&composition](auto in,auto out,auto b){composition.apply(in,out,b);},
                solver_options);
            if(prepare) {
                const auto rank=solver.rank();
                std::cout << "rank=" << rank.rank << " nullity=" << rank.nullity
                    << " handoff=" << solver_options.checkpoint_path.string()+".rank" << std::endl;
                solver.print_timing();return;
            }
            kernel=solver.nullspace();
            solver.print_timing();
            if(!kernel.complete)throw std::runtime_error("natural nullspace extraction incomplete");
            GraphHomology::convert_natural_representatives<fieldType>(kernel.basis,window.down,adjoints);
            saved_aut=std::move(adjoints.middle_aut);
        } else {
            const GraphHomology::StackedContraction stack(window.down, window.up);
            operator_rows=stack.rows();
            std::cout << "stack_transpose_bytes=" << stack.transpose_allocated_bytes()
                << " threads=" << stack.threads << std::endl;
            Solver solver(stack.rows(), stack.columns(),
                [&stack](auto in, auto out, auto b) { stack.apply<fieldType>(in,out,b); },
                [&stack](auto in, auto out, auto b) { stack.transpose<fieldType>(in,out,b); },
                solver_options);
            kernel=solver.nullspace();
            solver.print_timing();
            if (!kernel.complete) throw std::runtime_error("stack nullspace extraction incomplete");
            for (const auto& x : kernel.basis) {
                for(auto value:window.down.template apply<fieldType>(x))
                    if(value!=fieldType{})throw std::runtime_error("nonzero down residual");
                for (std::size_t c=0;c<window.up.columns();++c) {
                    fieldType value{};
                    for (auto i=window.up.offsets[c];i<window.up.offsets[c+1];++i)
                        value += x[window.up.row_indices[i]] * window.up.coefficients[i];
                    if (value != fieldType{}) throw std::runtime_error("nonzero up-transpose residual");
                }
            }
        }
        const std::filesystem::path destination(output);
        // A fresh folder prevents accidentally replacing saved representatives.
        if (!std::filesystem::create_directory(destination))
            throw std::runtime_error("nullspace output directory already exists");
        std::ofstream graphs(destination / "basis.tsv"), vectors(destination / "vectors.tsv"), metadata(destination / "metadata.txt");
        graphs << "basis_id\tordered_directed_edges\n";
        for (std::size_t i=0;i<window.middle.size();++i) {
            graphs << i;
            for (Int e=0;e<decltype(window.middle)::G::N_EDGES_;++e) {
                auto [a,b] = window.middle.graphs[i].getEdge(e);
                graphs << '\t' << int(a) << ',' << int(b);
            }
            graphs << '\n';
        }
        vectors << "vector_id\tbasis_id\tcoefficient_mod_p\n";
        for (std::size_t j=0;j<kernel.basis.size();++j)
            for (std::size_t i=0;i<kernel.basis[j].size();++i)
                if (kernel.basis[j][i]!=fieldType{}) vectors << j << '\t' << i << '\t' << kernel.basis[j][i].value() << '\n';
        metadata << "loop=" << L << " vertices=" << V << " parity=" << (P==GraphHomology::Parity::even?"even":"odd")
            << " prime=" << fieldType::characteristic()
            << (legacy ? "\noperator=[down;up^T] unweighted=true\n"
                : "\noperator=(S_down*C_down+C_up*S_up)^T internal_weights=false preconditioner=left_diagonal_after_reduction\nconversion=divide_by_aut_then_normalize residuals=C_down*x,S_up*x\n")
            << "rows=" << operator_rows << " columns=" << window.middle.size() << " nullity=" << kernel.basis.size()
            << "\ncomplete=true completeness_probabilistic=true residuals_verified=true\n"
            << "block_size=" << block_size << " seed=" << seed << " trials=1\n"
            << "orientation: vertices 0..V-1 in order; edges and their endpoint directions exactly as stored in basis.tsv\n";
        if(!legacy) {
            std::ofstream aut(destination / "automorphisms.tsv");
            aut << "basis_id\tautomorphism_order\n";
            for(std::size_t i=0;i<saved_aut.size();++i)aut << i << '\t' << saved_aut[i] << '\n';
            aut.close();
            if(!aut)throw std::runtime_error("failed writing automorphism orders");
        }
        graphs.close(); vectors.close(); metadata.close();
        if (!graphs || !vectors || !metadata) throw std::runtime_error("failed writing nullspace output");
        std::cout << "operator_rank_estimate=" << kernel.rank_estimate.rank << " nullspace_dimension=" << kernel.basis.size()
            << " residuals_verified=true nullspace_seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()
            << " output=" << output << std::endl;
    }
    if (ranks) {
        using Solver = VectorSpace::block_wiedemann_solver<fieldType>;
        const auto rank_start = std::chrono::steady_clock::now();
        Solver::options config;
        config.block_size = block_size;
        config.seed = seed;
        config.log = &std::cout;
        const auto down = Solver(window.down, config).rank();
        ++config.seed;
        const auto up = Solver(window.up, config).rank();
        const auto report = [](const char* name, const auto& result) {
            std::cout << name << " rank=" << result.rank << " nullity=" << result.nullity << " trials=";
            for (auto rank : result.trial_ranks) std::cout << rank << ',';
            std::cout << '\n';
        };
        report("down", down); report("up", up);
        if (down.nullity < up.rank) throw std::runtime_error("inconsistent rank estimates; retry");
        std::cout << "homology_dimension_estimate=" << down.nullity - up.rank
            << " method=block_wiedemann probabilistic=true block_size=" << block_size
            << " rank_seconds=" << std::chrono::duration<double>(
                std::chrono::steady_clock::now() - rank_start).count() << std::endl;
    }
}
int main(int argc, char** argv) {
    if (argc < 3 || argc > 7) { std::cerr << "Usage: " << argv[0] << " GRAPH_DIRECTORY even|odd [--rank [BLOCK_SIZE [SEED]] | --natural-adjoint [BLOCK_SIZE [SEED]] | --prepare-nullspace OUTPUT_DIRECTORY [BLOCK_SIZE [SEED]] | --nullspace OUTPUT_DIRECTORY [BLOCK_SIZE [SEED]] | --legacy-nullspace OUTPUT_DIRECTORY [BLOCK_SIZE [SEED]]]\n"; return 2; }
    try {
        const std::string parity(argv[2]);
        const bool ranks = argc >= 4 && std::string(argv[3]) == "--rank";
        const bool natural = argc >= 4 && std::string(argv[3]) == "--natural-adjoint";
        const bool legacy = argc >= 4 && std::string(argv[3]) == "--legacy-nullspace";
        const bool prepare = argc >= 4 && std::string(argv[3]) == "--prepare-nullspace";
        const bool nullspace = prepare || legacy || (argc >= 4 && std::string(argv[3]) == "--nullspace");
        if (argc >= 4 && !ranks && !nullspace && !natural) throw std::invalid_argument("expected --rank, --nullspace or --natural-adjoint");
        if (((ranks || natural) && argc > 6) || (nullspace && argc < 5)) throw std::invalid_argument("invalid option arguments");
        const char* output = nullspace ? argv[4] : nullptr;
        const int block_arg = nullspace ? 5 : 4;
        std::size_t block_size = 8;
        if (argc > block_arg) {
            const std::string value(argv[block_arg]);
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::invalid_argument("block size must be a positive integer");
            block_size = std::stoull(value);
            if (!block_size) throw std::invalid_argument("block size must be positive");
        }
        std::uint64_t seed = 17;
        if (argc > block_arg+1) {
            const std::string value(argv[block_arg+1]);
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::invalid_argument("seed must be an unsigned integer");
            seed = std::stoull(value);
        }
        if (ranks || nullspace)
            std::cout << "block_size=" << block_size << " seed=" << seed << " trials=1" << std::endl;
        if (parity == "even") run<GraphHomology::Parity::even>(argv[1], ranks, block_size, output, seed, natural, legacy, prepare);
        else if (parity == "odd") run<GraphHomology::Parity::odd>(argv[1], ranks, block_size, output, seed, natural, legacy, prepare);
        else throw std::invalid_argument("parity must be even or odd");
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
