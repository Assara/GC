# ======= CONFIGURATION =======

CXX       := g++
LD        := g++

TARGET    := gc
PLAYGROUND_TARGET := gc_contraction_playground
SPLIT_WATERFALL_TARGET := gc_split_waterfall_playground
STANDARDIZE_PERF_TARGET := standardize_perf
STANDARDIZE3_PROFILE_TARGET := standardize3_profile
ASSIGN_PERMUTED_SORT_COMPARE_TARGET := assign_permuted_sort_compare
GRAPH_STANDARDIZER_COMPARE_TARGET := graph_standardizer_compare
RANDOMIZE_SPLIT_CONTRACT_TARGET := randomize_split_contract_graphs
FILTER_CLASS_BY_VALENCE_TARGET := filter_class_by_valence
WHEEL_SPLIT_CONTRACT_REPS_TARGET := wheel_split_contract_reps
GRAPH_STAGE_GENERATOR_LOOP ?= 7
GRAPH_GENERATION_PROFILE ?= 0
GRAPH_GENERATION_PROFILE_SUFFIX := $(if $(filter 1,$(GRAPH_GENERATION_PROFILE)),_profile,)
GRAPH_STAGE_GENERATOR_TARGET := graph_stage_generator_$(GRAPH_STAGE_GENERATOR_LOOP)$(GRAPH_GENERATION_PROFILE_SUFFIX)
GENG_FINALIZER_TARGET := geng_finalizer_$(GRAPH_STAGE_GENERATOR_LOOP)
TRANSIENT_GRAPH_TEST_TARGET := transient_graph_test
ROOTED_TRANSIENT_TEST_TARGET := rooted_transient_test
AGGREGATE_TRANSIENT_TEST_TARGET := aggregate_transient_test
SUPPORT_TRANSIENT_TEST_TARGET := support_transient_test
MAPPED_SUPPORT_TRANSIENT_TEST_TARGET := mapped_support_transient_test
FINAL_CANONICALIZATION_TEST_TARGET := final_canonicalization_test
LINEAR_PROBE_SET_TEST_TARGET := linear_probe_set_test
TRICONNECTED_ORACLE_LOOP ?= 8
TRICONNECTED_ORACLE_TARGET := compare_triconnected_generation_$(TRICONNECTED_ORACLE_LOOP)
MIN_TADPOLES ?= 3
MAX_TADPOLES ?= 5
DIMENSION_DIR ?= output/graph_dimensions
DIMENSION_TABLE ?= $(DIMENSION_DIR)/dimensions.tsv
SOLVER_COMPARE_TARGET := solver_compare
TEST_TARGET := gc_test
BUILD_DIR := build

WHEEL ?= 9
SPLIT_CONTRACT_ROUNDS ?= $(shell echo $$((($(WHEEL) - 3) / 2)))
SPLIT_CONTRACT_DIR ?= output/split_contract
SPLIT_CONTRACT_CANDIDATES ?= $(SPLIT_CONTRACT_DIR)/W$(WHEEL)_candidates.txt
SPLIT_CONTRACT_MAP_PREFIX ?= $(SPLIT_CONTRACT_DIR)/maps/W$(WHEEL)_rounds$(SPLIT_CONTRACT_ROUNDS)_split_map
CHECKPOINT_INTERVAL ?= 32

# Include directories
INC       := -I. -IVectorSpace

# ---- OpenMP ----
# GCC uses the installed libgomp OpenMP runtime.
OPENMP_CXXFLAGS := -fopenmp
OPENMP_LDFLAGS  := -fopenmp

# Prefer LLD when installed, otherwise let the compiler select the system linker.
LINKER_FLAG := $(if $(shell command -v ld.lld 2>/dev/null),-fuse-ld=lld,)

NAUTY_CXXFLAGS ?= $(shell pkg-config --cflags nauty 2>/dev/null)
NAUTY_LDLIBS ?= $(shell pkg-config --libs nauty 2>/dev/null || echo -lnauty)
LINBOX_CXXFLAGS ?= $(shell pkg-config --cflags linbox givaro fflas-ffpack 2>/dev/null)
LINBOX_LDLIBS ?= $(shell pkg-config --libs linbox givaro fflas-ffpack 2>/dev/null || echo -llinbox -lgivaro -lfflas -lffpack -lgmp)

# Compiler and linker flags
CXXFLAGS  := -std=c++23 -O3 -march=native -Wall -Wextra -Wpedantic $(INC) $(OPENMP_CXXFLAGS)
LDFLAGS   := $(LINKER_FLAG) $(OPENMP_LDFLAGS)

# ======= SOURCE / BUILD SETUP =======

# Main binary sources
MAIN_SRCS := main.cpp GC_split_playground.cpp
MAIN_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(MAIN_SRCS))

# Standalone binaries
PLAYGROUND_SRCS := GC_contraction_playground.cpp
PLAYGROUND_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(PLAYGROUND_SRCS))

SPLIT_WATERFALL_SRCS := GC_split_waterfall_playground.cpp
SPLIT_WATERFALL_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SPLIT_WATERFALL_SRCS))

STANDARDIZE_PERF_SRCS := tools/standardize_perf.cpp
STANDARDIZE_PERF_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(STANDARDIZE_PERF_SRCS))

STANDARDIZE3_PROFILE_SRCS := tools/standardize3_profile.cpp
STANDARDIZE3_PROFILE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(STANDARDIZE3_PROFILE_SRCS))

ASSIGN_PERMUTED_SORT_COMPARE_SRCS := tools/assign_permuted_sort_compare.cpp
ASSIGN_PERMUTED_SORT_COMPARE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(ASSIGN_PERMUTED_SORT_COMPARE_SRCS))

GRAPH_STANDARDIZER_COMPARE_SRCS := tools/graph_standardizer_comparison.cpp
GRAPH_STANDARDIZER_COMPARE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(GRAPH_STANDARDIZER_COMPARE_SRCS))
GRAPH_STANDARDIZER_COMPARE_NAUTY_OBJS := $(BUILD_DIR)/tools/graph_standardizer_comparison_nauty.o

RANDOMIZE_SPLIT_CONTRACT_SRCS := tools/randomize_split_contract_graphs.cpp
RANDOMIZE_SPLIT_CONTRACT_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(RANDOMIZE_SPLIT_CONTRACT_SRCS))

FILTER_CLASS_BY_VALENCE_SRCS := tools/filter_class_by_valence.cpp
FILTER_CLASS_BY_VALENCE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(FILTER_CLASS_BY_VALENCE_SRCS))

WHEEL_SPLIT_CONTRACT_REPS_SRCS := tools/wheel_split_contract_reps.cpp
WHEEL_SPLIT_CONTRACT_REPS_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(WHEEL_SPLIT_CONTRACT_REPS_SRCS))

GRAPH_STAGE_GENERATOR_SRC := tools/generate_graph_stages.cpp
GRAPH_STAGE_GENERATOR_OBJ := $(BUILD_DIR)/tools/generate_graph_stages_L$(GRAPH_STAGE_GENERATOR_LOOP)$(GRAPH_GENERATION_PROFILE_SUFFIX).o
GENG_FINALIZER_SRC := tools/finalize_geng_stream.cpp
GENG_FINALIZER_OBJ := $(BUILD_DIR)/tools/finalize_geng_stream_L$(GRAPH_STAGE_GENERATOR_LOOP).o
TRANSIENT_GRAPH_TEST_SRC := tools/transient_graph_test.cpp
TRANSIENT_GRAPH_TEST_OBJ := $(BUILD_DIR)/tools/transient_graph_test.o
ROOTED_TRANSIENT_TEST_SRC := tools/rooted_transient_test.cpp
ROOTED_TRANSIENT_TEST_OBJ := $(BUILD_DIR)/tools/rooted_transient_test.o
AGGREGATE_TRANSIENT_TEST_SRC := tools/aggregate_transient_test.cpp
AGGREGATE_TRANSIENT_TEST_OBJ := $(BUILD_DIR)/tools/aggregate_transient_test.o
SUPPORT_TRANSIENT_TEST_SRC := tools/support_transient_test.cpp
SUPPORT_TRANSIENT_TEST_OBJ := $(BUILD_DIR)/tools/support_transient_test.o
MAPPED_SUPPORT_TRANSIENT_TEST_SRC := tools/mapped_support_transient_test.cpp
MAPPED_SUPPORT_TRANSIENT_TEST_OBJ := $(BUILD_DIR)/tools/mapped_support_transient_test.o
FINAL_CANONICALIZATION_TEST_SRC := tools/final_canonicalization_test.cpp
FINAL_CANONICALIZATION_TEST_OBJ := $(BUILD_DIR)/tools/final_canonicalization_test.o
LINEAR_PROBE_SET_TEST_SRC := tools/linear_probe_set_test.cpp
LINEAR_PROBE_SET_TEST_OBJ := $(BUILD_DIR)/tools/linear_probe_set_test.o
TRICONNECTED_ORACLE_SRC := tools/compare_triconnected_generation.cpp
TRICONNECTED_ORACLE_OBJ := $(BUILD_DIR)/tools/compare_triconnected_generation_L$(TRICONNECTED_ORACLE_LOOP).o
GRAPH_GENERATION_HEADERS := GraphGeneration/SupportTransientGraph.hpp GraphGeneration/SupportTransientStandardizer.hpp GraphGeneration/UnrootedSupportTransientGraph.hpp GraphGeneration/UnrootedSupportTransientStandardizer.hpp GraphGeneration/RootedTransientGraph.hpp GraphGeneration/MappedSupportTransientFile.hpp GraphGeneration/FinalCanonicalization.hpp GraphGeneration/MappedGraphFile.hpp GraphGeneration/MappedFinalGraphFile.hpp CombinatorialUtils.hpp LinearProbeSet.hpp graph.hpp GraphStandardizer.hpp

SOLVER_COMPARE_SRCS := tools/solver_comparison.cpp
SOLVER_COMPARE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOLVER_COMPARE_SRCS))
SOLVER_COMPARE_LINBOX_OBJS := $(BUILD_DIR)/tools/solver_comparison_linbox.o

TEST_SRCS := GC_test.cpp
TEST_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SRCS))

# ======= RULES =======

.PHONY: all clean run playground split-waterfall run-split-waterfall standardize-perf standardize3-profile compare-standardize3-commits assign-permuted-sort-compare randomize-split-contract-graphs filter-class-by-valence graph-stage-generator geng-finalizer graph-dimension-table transient-graph-test run-transient-graph-test rooted-transient-test run-rooted-transient-test aggregate-transient-test run-aggregate-transient-test support-transient-test run-support-transient-test mapped-support-transient-test run-mapped-support-transient-test final-canonicalization-test run-final-canonicalization-test linear-probe-set-test run-linear-probe-set-test triconnected-oracle run-sparse-rank-test run-standardize-perf run-standardize3-profile run-compare-standardize3-commits run-assign-permuted-sort-compare run-randomize-split-contract-graphs run-filter-class-by-valence graph-standardizer-compare graph-standardizer-compare-nauty run-graph-standardizer-compare wheel-split-contract-reps run-wheel-split-contract-reps run-wheel-class run-wheel-class-generated-constrained run-split-contract-files run-split-contract-solver solver-compare solver-compare-linbox run-solver-compare test

all: $(TARGET)

$(TARGET): $(MAIN_OBJS)
	@echo "🚧 Linking $(TARGET) with GCC..."
	$(LD) $(MAIN_OBJS) $(LDFLAGS) -o $(TARGET)
	@echo "✅ Build complete."

$(PLAYGROUND_TARGET): $(PLAYGROUND_OBJS)
	@echo "🚧 Linking $(PLAYGROUND_TARGET) with Clang + LLD..."
	$(LD) $(PLAYGROUND_OBJS) $(LDFLAGS) -o $(PLAYGROUND_TARGET)
	@echo "✅ Build complete."

$(SPLIT_WATERFALL_TARGET): $(SPLIT_WATERFALL_OBJS)
	@echo "🚧 Linking $(SPLIT_WATERFALL_TARGET) with Clang + LLD..."
	$(LD) $(SPLIT_WATERFALL_OBJS) $(LDFLAGS) -o $(SPLIT_WATERFALL_TARGET)
	@echo "✅ Build complete."

$(STANDARDIZE_PERF_TARGET): $(STANDARDIZE_PERF_OBJS)
	@echo "🚧 Linking $(STANDARDIZE_PERF_TARGET) with Clang + LLD..."
	$(LD) $(STANDARDIZE_PERF_OBJS) $(LDFLAGS) -o $(STANDARDIZE_PERF_TARGET)
	@echo "✅ Build complete."

$(STANDARDIZE3_PROFILE_TARGET): $(STANDARDIZE3_PROFILE_OBJS)
	@echo "🚧 Linking $(STANDARDIZE3_PROFILE_TARGET) with Clang + LLD..."
	$(LD) $(STANDARDIZE3_PROFILE_OBJS) $(LDFLAGS) -o $(STANDARDIZE3_PROFILE_TARGET)
	@echo "✅ Build complete."

$(ASSIGN_PERMUTED_SORT_COMPARE_TARGET): $(ASSIGN_PERMUTED_SORT_COMPARE_OBJS)
	@echo "🚧 Linking $(ASSIGN_PERMUTED_SORT_COMPARE_TARGET) with Clang + LLD..."
	$(LD) $(ASSIGN_PERMUTED_SORT_COMPARE_OBJS) $(LDFLAGS) -o $(ASSIGN_PERMUTED_SORT_COMPARE_TARGET)
	@echo "✅ Build complete."

$(GRAPH_STANDARDIZER_COMPARE_TARGET): $(GRAPH_STANDARDIZER_COMPARE_OBJS)
	@echo "🚧 Linking $(GRAPH_STANDARDIZER_COMPARE_TARGET) with Clang + LLD..."
	$(LD) $(GRAPH_STANDARDIZER_COMPARE_OBJS) $(LDFLAGS) -o $(GRAPH_STANDARDIZER_COMPARE_TARGET)
	@echo "✅ Build complete."

$(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty: $(GRAPH_STANDARDIZER_COMPARE_NAUTY_OBJS)
	@echo "🚧 Linking $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty with Clang + LLD..."
	$(LD) $(GRAPH_STANDARDIZER_COMPARE_NAUTY_OBJS) $(LDFLAGS) $(NAUTY_LDLIBS) -o $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty
	@echo "✅ Build complete."

$(RANDOMIZE_SPLIT_CONTRACT_TARGET): $(RANDOMIZE_SPLIT_CONTRACT_OBJS)
	@echo "🚧 Linking $(RANDOMIZE_SPLIT_CONTRACT_TARGET) with Clang + LLD..."
	$(LD) $(RANDOMIZE_SPLIT_CONTRACT_OBJS) $(LDFLAGS) -o $(RANDOMIZE_SPLIT_CONTRACT_TARGET)
	@echo "✅ Build complete."

$(FILTER_CLASS_BY_VALENCE_TARGET): $(FILTER_CLASS_BY_VALENCE_OBJS)
	@echo "🚧 Linking $(FILTER_CLASS_BY_VALENCE_TARGET) with Clang + LLD..."
	$(LD) $(FILTER_CLASS_BY_VALENCE_OBJS) $(LDFLAGS) -o $(FILTER_CLASS_BY_VALENCE_TARGET)
	@echo "✅ Build complete."

$(WHEEL_SPLIT_CONTRACT_REPS_TARGET): $(WHEEL_SPLIT_CONTRACT_REPS_OBJS)
	@echo "🚧 Linking $(WHEEL_SPLIT_CONTRACT_REPS_TARGET) with Clang + LLD..."
	$(LD) $(WHEEL_SPLIT_CONTRACT_REPS_OBJS) $(LDFLAGS) -o $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "✅ Build complete."

$(GRAPH_STAGE_GENERATOR_TARGET): $(GRAPH_STAGE_GENERATOR_OBJ)
	@echo "🚧 Linking $(GRAPH_STAGE_GENERATOR_TARGET) with Clang + LLD..."
	$(LD) $(GRAPH_STAGE_GENERATOR_OBJ) $(LDFLAGS) -o $(GRAPH_STAGE_GENERATOR_TARGET)
	@echo "✅ Build complete."

$(GRAPH_STAGE_GENERATOR_OBJ): $(GRAPH_STAGE_GENERATOR_SRC) $(GRAPH_GENERATION_HEADERS)
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $< for loop number $(GRAPH_STAGE_GENERATOR_LOOP)..."
	$(CXX) $(CXXFLAGS) -DGC_GENERATION_LOOP=$(GRAPH_STAGE_GENERATOR_LOOP) $(if $(filter 1,$(GRAPH_GENERATION_PROFILE)),-DGC_PROFILE_GRAPH_GENERATION) -c $< -o $@

$(GENG_FINALIZER_TARGET): $(GENG_FINALIZER_OBJ)
	@echo "🚧 Linking $(GENG_FINALIZER_TARGET) with Clang + LLD..."
	$(LD) $(GENG_FINALIZER_OBJ) $(LDFLAGS) -o $(GENG_FINALIZER_TARGET)
	@echo "✅ Build complete."

$(GENG_FINALIZER_OBJ): $(GENG_FINALIZER_SRC) $(GRAPH_GENERATION_HEADERS)
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $< for loop number $(GRAPH_STAGE_GENERATOR_LOOP)..."
	$(CXX) $(CXXFLAGS) -DGC_GENERATION_LOOP=$(GRAPH_STAGE_GENERATOR_LOOP) -c $< -o $@

$(TRANSIENT_GRAPH_TEST_TARGET): $(TRANSIENT_GRAPH_TEST_OBJ)
	@echo "🚧 Linking $(TRANSIENT_GRAPH_TEST_TARGET) with Clang + LLD..."
	$(LD) $(TRANSIENT_GRAPH_TEST_OBJ) $(LDFLAGS) -o $(TRANSIENT_GRAPH_TEST_TARGET)
	@echo "✅ Build complete."

$(TRANSIENT_GRAPH_TEST_OBJ): $(TRANSIENT_GRAPH_TEST_SRC) GraphGeneration/TransientGraph.hpp CombinatorialUtils.hpp

$(ROOTED_TRANSIENT_TEST_TARGET): $(ROOTED_TRANSIENT_TEST_OBJ)
	@echo "🚧 Linking $(ROOTED_TRANSIENT_TEST_TARGET) with Clang + LLD..."
	$(LD) $(ROOTED_TRANSIENT_TEST_OBJ) $(LDFLAGS) -o $(ROOTED_TRANSIENT_TEST_TARGET)
	@echo "✅ Build complete."

$(ROOTED_TRANSIENT_TEST_OBJ): $(ROOTED_TRANSIENT_TEST_SRC) GraphGeneration/RootedTransientGraph.hpp GraphGeneration/RootedTransientStandardizer.hpp GraphGeneration/TransientGraph.hpp CombinatorialUtils.hpp

$(AGGREGATE_TRANSIENT_TEST_TARGET): $(AGGREGATE_TRANSIENT_TEST_OBJ)
	@echo "🚧 Linking $(AGGREGATE_TRANSIENT_TEST_TARGET) with Clang + LLD..."
	$(LD) $(AGGREGATE_TRANSIENT_TEST_OBJ) $(LDFLAGS) -o $(AGGREGATE_TRANSIENT_TEST_TARGET)
	@echo "✅ Build complete."

$(AGGREGATE_TRANSIENT_TEST_OBJ): $(AGGREGATE_TRANSIENT_TEST_SRC) GraphGeneration/AggregateTransientGraph.hpp GraphGeneration/AggregateTransientStandardizer.hpp GraphGeneration/RootedTransientGraph.hpp GraphGeneration/RootedTransientStandardizer.hpp GraphGeneration/TransientGraph.hpp CombinatorialUtils.hpp

$(SUPPORT_TRANSIENT_TEST_TARGET): $(SUPPORT_TRANSIENT_TEST_OBJ)
	@echo "🚧 Linking $(SUPPORT_TRANSIENT_TEST_TARGET) with Clang + LLD..."
	$(LD) $(SUPPORT_TRANSIENT_TEST_OBJ) $(LDFLAGS) -o $(SUPPORT_TRANSIENT_TEST_TARGET)
	@echo "✅ Build complete."

$(SUPPORT_TRANSIENT_TEST_OBJ): $(SUPPORT_TRANSIENT_TEST_SRC) GraphGeneration/SupportTransientGraph.hpp GraphGeneration/SupportTransientStandardizer.hpp GraphGeneration/RootedTransientGraph.hpp GraphGeneration/RootedTransientStandardizer.hpp GraphGeneration/TransientGraph.hpp CombinatorialUtils.hpp

$(MAPPED_SUPPORT_TRANSIENT_TEST_TARGET): $(MAPPED_SUPPORT_TRANSIENT_TEST_OBJ)
	@echo "🚧 Linking $(MAPPED_SUPPORT_TRANSIENT_TEST_TARGET) with Clang + LLD..."
	$(LD) $(MAPPED_SUPPORT_TRANSIENT_TEST_OBJ) $(LDFLAGS) -o $(MAPPED_SUPPORT_TRANSIENT_TEST_TARGET)
	@echo "✅ Build complete."

$(MAPPED_SUPPORT_TRANSIENT_TEST_OBJ): $(MAPPED_SUPPORT_TRANSIENT_TEST_SRC) GraphGeneration/MappedSupportTransientFile.hpp GraphGeneration/MappedGraphFile.hpp GraphGeneration/SupportTransientGraph.hpp

$(FINAL_CANONICALIZATION_TEST_TARGET): $(FINAL_CANONICALIZATION_TEST_OBJ)
	@echo "🚧 Linking $(FINAL_CANONICALIZATION_TEST_TARGET) with Clang + LLD..."
	$(LD) $(FINAL_CANONICALIZATION_TEST_OBJ) $(LDFLAGS) -o $(FINAL_CANONICALIZATION_TEST_TARGET)
	@echo "✅ Build complete."

$(FINAL_CANONICALIZATION_TEST_OBJ): $(FINAL_CANONICALIZATION_TEST_SRC) GraphGeneration/FinalCanonicalization.hpp

$(LINEAR_PROBE_SET_TEST_TARGET): $(LINEAR_PROBE_SET_TEST_OBJ)
	@echo "🚧 Linking $(LINEAR_PROBE_SET_TEST_TARGET) with Clang + LLD..."
	$(LD) $(LINEAR_PROBE_SET_TEST_OBJ) $(LDFLAGS) -o $(LINEAR_PROBE_SET_TEST_TARGET)
	@echo "✅ Build complete."

$(LINEAR_PROBE_SET_TEST_OBJ): $(LINEAR_PROBE_SET_TEST_SRC) LinearProbeSet.hpp GraphGeneration/SupportTransientGraph.hpp graph.hpp graph_hash.hpp

$(TRICONNECTED_ORACLE_TARGET): $(TRICONNECTED_ORACLE_OBJ)
	@echo "🚧 Linking $(TRICONNECTED_ORACLE_TARGET) with Clang + LLD..."
	$(LD) $(TRICONNECTED_ORACLE_OBJ) $(LDFLAGS) -o $(TRICONNECTED_ORACLE_TARGET)
	@echo "✅ Build complete."

$(TRICONNECTED_ORACLE_OBJ): $(TRICONNECTED_ORACLE_SRC) GraphGeneration/FinalCanonicalization.hpp GraphGeneration/MappedFinalGraphFile.hpp GraphGeneration/MappedGraphFile.hpp graph.hpp GraphStandardizer.hpp
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $< for loop number $(TRICONNECTED_ORACLE_LOOP)..."
	$(CXX) $(CXXFLAGS) -DGC_TRICONNECTED_ORACLE_LOOP=$(TRICONNECTED_ORACLE_LOOP) -c $< -o $@

$(SOLVER_COMPARE_TARGET): $(SOLVER_COMPARE_OBJS)
	@echo "🚧 Linking $(SOLVER_COMPARE_TARGET) with Clang + LLD..."
	$(LD) $(SOLVER_COMPARE_OBJS) $(LDFLAGS) -o $(SOLVER_COMPARE_TARGET)
	@echo "✅ Build complete."

$(SOLVER_COMPARE_TARGET)-linbox: $(SOLVER_COMPARE_LINBOX_OBJS)
	@echo "🚧 Linking $(SOLVER_COMPARE_TARGET)-linbox with Clang + LLD..."
	$(LD) $(SOLVER_COMPARE_LINBOX_OBJS) $(LDFLAGS) $(LINBOX_LDLIBS) -o $(SOLVER_COMPARE_TARGET)-linbox
	@echo "✅ Build complete."

$(TEST_TARGET): $(TEST_OBJS)
	@echo "🚧 Linking $(TEST_TARGET) with Clang + LLD..."
	$(LD) $(TEST_OBJS) $(LDFLAGS) -o $(TEST_TARGET)
	@echo "✅ Build complete."

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $<..."
	$(CXX) $(CXXFLAGS) $(if $(filter tools/graph_standardizer_comparison.cpp tools/standardize3_profile.cpp,$<),-DGC_PROFILE_STANDARDIZER_SORT) -c $< -o $@

$(GRAPH_STANDARDIZER_COMPARE_NAUTY_OBJS): tools/graph_standardizer_comparison.cpp
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $< with Nauty support..."
	$(CXX) $(CXXFLAGS) $(NAUTY_CXXFLAGS) -DGC_PERF_WITH_NAUTY -DGC_PROFILE_STANDARDIZER_SORT -c $< -o $@

$(SOLVER_COMPARE_LINBOX_OBJS): tools/solver_comparison.cpp
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $< with LinBox support..."
	$(CXX) $(CXXFLAGS) $(LINBOX_CXXFLAGS) -DGC_PERF_WITH_LINBOX -c $< -o $@

clean:
	@echo "🧹 Cleaning build files..."
	rm -rf $(BUILD_DIR) $(TARGET) $(PLAYGROUND_TARGET) $(SPLIT_WATERFALL_TARGET) $(STANDARDIZE_PERF_TARGET) $(STANDARDIZE3_PROFILE_TARGET) $(ASSIGN_PERMUTED_SORT_COMPARE_TARGET) $(GRAPH_STANDARDIZER_COMPARE_TARGET) $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty $(RANDOMIZE_SPLIT_CONTRACT_TARGET) $(FILTER_CLASS_BY_VALENCE_TARGET) $(WHEEL_SPLIT_CONTRACT_REPS_TARGET) graph_stage_generator_* geng_finalizer_* $(TRANSIENT_GRAPH_TEST_TARGET) $(ROOTED_TRANSIENT_TEST_TARGET) $(AGGREGATE_TRANSIENT_TEST_TARGET) $(SUPPORT_TRANSIENT_TEST_TARGET) $(MAPPED_SUPPORT_TRANSIENT_TEST_TARGET) $(FINAL_CANONICALIZATION_TEST_TARGET) $(LINEAR_PROBE_SET_TEST_TARGET) compare_triconnected_generation_* $(SOLVER_COMPARE_TARGET) $(SOLVER_COMPARE_TARGET)-linbox $(TEST_TARGET)

run: all
	@./$(TARGET)

playground: $(PLAYGROUND_TARGET)

split-waterfall: $(SPLIT_WATERFALL_TARGET)

run-split-waterfall: $(SPLIT_WATERFALL_TARGET)
	@./$(SPLIT_WATERFALL_TARGET) $(or $(WHEEL),7) $(or $(OUT),W$(or $(WHEEL),7).txt)

standardize-perf: $(STANDARDIZE_PERF_TARGET)

standardize3-profile: $(STANDARDIZE3_PROFILE_TARGET)

compare-standardize3-commits:

assign-permuted-sort-compare: $(ASSIGN_PERMUTED_SORT_COMPARE_TARGET)

randomize-split-contract-graphs: $(RANDOMIZE_SPLIT_CONTRACT_TARGET)

filter-class-by-valence: $(FILTER_CLASS_BY_VALENCE_TARGET)

graph-stage-generator: $(GRAPH_STAGE_GENERATOR_TARGET)

geng-finalizer: $(GENG_FINALIZER_TARGET)

graph-dimension-table:
	@MIN_TADPOLES="$(MIN_TADPOLES)" ./tools/generate_dimension_table.sh "$(MAX_TADPOLES)" "$(DIMENSION_DIR)" "$(DIMENSION_TABLE)"

transient-graph-test: $(TRANSIENT_GRAPH_TEST_TARGET)

.PHONY: transient-graph2-test run-transient-graph2-test
transient-graph2-test: $(BUILD_DIR)/transient_graph2_test

$(BUILD_DIR)/transient_graph2_test: $(BUILD_DIR)/tools/transient_graph2_test.o
	$(LD) $< $(LDFLAGS) -o $@

$(BUILD_DIR)/tools/transient_graph2_test.o: tools/transient_graph2_test.cpp GraphGeneration/TransientGraph2.hpp GraphGeneration/TransientGraph2Standardizer.hpp graph.hpp GraphStandardizer.hpp

run-transient-graph2-test: $(BUILD_DIR)/transient_graph2_test
	@./$(BUILD_DIR)/transient_graph2_test

run-transient-graph-test: $(TRANSIENT_GRAPH_TEST_TARGET)
	@./$(TRANSIENT_GRAPH_TEST_TARGET)

rooted-transient-test: $(ROOTED_TRANSIENT_TEST_TARGET)

run-rooted-transient-test: $(ROOTED_TRANSIENT_TEST_TARGET)
	@./$(ROOTED_TRANSIENT_TEST_TARGET)

aggregate-transient-test: $(AGGREGATE_TRANSIENT_TEST_TARGET)

run-aggregate-transient-test: $(AGGREGATE_TRANSIENT_TEST_TARGET)
	@./$(AGGREGATE_TRANSIENT_TEST_TARGET)

support-transient-test: $(SUPPORT_TRANSIENT_TEST_TARGET)

run-support-transient-test: $(SUPPORT_TRANSIENT_TEST_TARGET)
	@./$(SUPPORT_TRANSIENT_TEST_TARGET)

mapped-support-transient-test: $(MAPPED_SUPPORT_TRANSIENT_TEST_TARGET)

run-mapped-support-transient-test: $(MAPPED_SUPPORT_TRANSIENT_TEST_TARGET)
	@./$(MAPPED_SUPPORT_TRANSIENT_TEST_TARGET)

final-canonicalization-test: $(FINAL_CANONICALIZATION_TEST_TARGET)

run-final-canonicalization-test: $(FINAL_CANONICALIZATION_TEST_TARGET)
	@./$(FINAL_CANONICALIZATION_TEST_TARGET)

linear-probe-set-test: $(LINEAR_PROBE_SET_TEST_TARGET)

run-linear-probe-set-test: $(LINEAR_PROBE_SET_TEST_TARGET)
	@./$(LINEAR_PROBE_SET_TEST_TARGET)

triconnected-oracle: $(TRICONNECTED_ORACLE_TARGET)

run-sparse-rank-test:
	@$(MAKE) -C VectorSpace/tests "$(CURDIR)/VectorSpace/tests/bin/test_sparse_rank"
	@./VectorSpace/tests/bin/test_sparse_rank

run-standardize-perf: $(STANDARDIZE_PERF_TARGET)
	@./$(STANDARDIZE_PERF_TARGET) $(or $(WHEEL),33) $(or $(REPEAT),1) $(or $(ITER),3)

run-standardize3-profile: $(STANDARDIZE3_PROFILE_TARGET)
	@./$(STANDARDIZE3_PROFILE_TARGET) $(or $(WHEEL),11) $(or $(ROUNDS),2) $(or $(REPEAT),1) $(or $(ITER),3) $(or $(WORKLOAD),split-contract)

run-compare-standardize3-commits:
	@./tools/compare_standardize3_commits.sh $(or $(OLD),d11ef3d) $(or $(NEW),HEAD) $(or $(WHEEL),11) $(or $(ROUNDS),1) $(or $(REPEAT),1) $(or $(ITER),3) $(or $(WORKLOAD),split-contract) $(FILE)

run-assign-permuted-sort-compare: $(ASSIGN_PERMUTED_SORT_COMPARE_TARGET)
	@./$(ASSIGN_PERMUTED_SORT_COMPARE_TARGET) $(or $(WHEEL),25) $(or $(ROUNDS),0) $(or $(ITER),3) $(or $(WORKLOAD),vmax)

run-randomize-split-contract-graphs: $(RANDOMIZE_SPLIT_CONTRACT_TARGET)
	@./$(RANDOMIZE_SPLIT_CONTRACT_TARGET) $(or $(WHEEL),11) $(IN) $(OUT) $(or $(SEED),123456789)

run-filter-class-by-valence: $(FILTER_CLASS_BY_VALENCE_TARGET)
	@./$(FILTER_CLASS_BY_VALENCE_TARGET) $(or $(WHEEL),7) $(IN) $(OUT)

graph-standardizer-compare: $(GRAPH_STANDARDIZER_COMPARE_TARGET)

graph-standardizer-compare-nauty: $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty

run-graph-standardizer-compare: $(GRAPH_STANDARDIZER_COMPARE_TARGET)
	@./$(GRAPH_STANDARDIZER_COMPARE_TARGET) $(or $(WHEEL),25) $(or $(ROUNDS),2) $(or $(REPEAT),1) $(or $(ITER),3) $(or $(WORKLOAD),split-contract)

wheel-split-contract-reps: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)

run-wheel-split-contract-reps: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(or $(WHEEL),11) $(or $(ROUNDS),2) $(OUT)

run-wheel-class: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "Generating wheel class via same-degree split-contract support + Wiedemann solve"
	@echo "Representatives output: $(OUT)"
	@echo "Class output: $(CLASS_OUT)"
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) $(OUT) $(CLASS_OUT)

run-wheel-class-generated-constrained: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "Generating wheel class via generated constrained split-contract solve"
	@echo "Representatives output: $(OUT)"
	@echo "Class output: $(CLASS_OUT)"
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) $(OUT) $(CLASS_OUT) generated-constrained

run-split-contract-files: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "Generating split-contract files for W$(WHEEL), rounds $(SPLIT_CONTRACT_ROUNDS)"
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) $(SPLIT_CONTRACT_CANDIDATES) unused enumerate
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) $(SPLIT_CONTRACT_CANDIDATES) $(SPLIT_CONTRACT_MAP_PREFIX) split-map

run-split-contract-solver: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "Solving from generated split-contract files for W$(WHEEL), rounds $(SPLIT_CONTRACT_ROUNDS)"
	@echo "This writes solver coefficients only; it does not reconstruct/save the final class."
	@echo "Checkpoint prefix: $(SPLIT_CONTRACT_MAP_PREFIX)_wiedemann_checkpoint.bin"
	@echo "Checkpoint interval: $(CHECKPOINT_INTERVAL) MMT iterations"
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) unused $(SPLIT_CONTRACT_MAP_PREFIX) split-map-solve $(CHECKPOINT_INTERVAL)

solver-compare: $(SOLVER_COMPARE_TARGET)

solver-compare-linbox: $(SOLVER_COMPARE_TARGET)-linbox

run-solver-compare: $(SOLVER_COMPARE_TARGET)
	@./$(SOLVER_COMPARE_TARGET) $(or $(ROWS),5000) $(or $(COLS),5000) $(or $(NNZ),50000) $(or $(SEED),17)

test: $(TEST_TARGET)
