# ======= CONFIGURATION =======

CXX       := clang++
LD        := clang++

TARGET    := gc
PLAYGROUND_TARGET := gc_contraction_playground
SPLIT_WATERFALL_TARGET := gc_split_waterfall_playground
STANDARDIZE_PERF_TARGET := standardize_perf
STANDARDIZE3_PROFILE_TARGET := standardize3_profile
ASSIGN_PERMUTED_SORT_COMPARE_TARGET := assign_permuted_sort_compare
GRAPH_STANDARDIZER_COMPARE_TARGET := graph_standardizer_compare
RANDOMIZE_SPLIT_CONTRACT_TARGET := randomize_split_contract_graphs
WHEEL_SPLIT_CONTRACT_REPS_TARGET := wheel_split_contract_reps
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
# For Clang on Linux, libomp is typical. If your system uses libgomp instead, see notes below.
OPENMP_CXXFLAGS := -fopenmp=libomp
OPENMP_LDFLAGS  := -fopenmp=libomp -lomp

NAUTY_CXXFLAGS ?= $(shell pkg-config --cflags nauty 2>/dev/null)
NAUTY_LDLIBS ?= $(shell pkg-config --libs nauty 2>/dev/null || echo -lnauty)
LINBOX_CXXFLAGS ?= $(shell pkg-config --cflags linbox givaro fflas-ffpack 2>/dev/null)
LINBOX_LDLIBS ?= $(shell pkg-config --libs linbox givaro fflas-ffpack 2>/dev/null || echo -llinbox -lgivaro -lfflas -lffpack -lgmp)

# Compiler and linker flags
CXXFLAGS  := -std=c++23 -O3 -march=native -Wall -Wextra -Wpedantic $(INC) -flto $(OPENMP_CXXFLAGS)
LDFLAGS   := -flto -fuse-ld=lld $(OPENMP_LDFLAGS)

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

WHEEL_SPLIT_CONTRACT_REPS_SRCS := tools/wheel_split_contract_reps.cpp
WHEEL_SPLIT_CONTRACT_REPS_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(WHEEL_SPLIT_CONTRACT_REPS_SRCS))

SOLVER_COMPARE_SRCS := tools/solver_comparison.cpp
SOLVER_COMPARE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOLVER_COMPARE_SRCS))
SOLVER_COMPARE_LINBOX_OBJS := $(BUILD_DIR)/tools/solver_comparison_linbox.o

TEST_SRCS := GC_test.cpp
TEST_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SRCS))

# ======= RULES =======

.PHONY: all clean run playground split-waterfall run-split-waterfall standardize-perf standardize3-profile compare-standardize3-commits assign-permuted-sort-compare randomize-split-contract-graphs run-standardize-perf run-standardize3-profile run-compare-standardize3-commits run-assign-permuted-sort-compare run-randomize-split-contract-graphs graph-standardizer-compare graph-standardizer-compare-nauty run-graph-standardizer-compare wheel-split-contract-reps run-wheel-split-contract-reps run-split-contract-files run-split-contract-solver solver-compare solver-compare-linbox run-solver-compare test

all: $(TARGET)

$(TARGET): $(MAIN_OBJS)
	@echo "🚧 Linking $(TARGET) with Clang + LLD..."
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

$(WHEEL_SPLIT_CONTRACT_REPS_TARGET): $(WHEEL_SPLIT_CONTRACT_REPS_OBJS)
	@echo "🚧 Linking $(WHEEL_SPLIT_CONTRACT_REPS_TARGET) with Clang + LLD..."
	$(LD) $(WHEEL_SPLIT_CONTRACT_REPS_OBJS) $(LDFLAGS) -o $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "✅ Build complete."

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
	rm -rf $(BUILD_DIR) $(TARGET) $(PLAYGROUND_TARGET) $(SPLIT_WATERFALL_TARGET) $(STANDARDIZE_PERF_TARGET) $(STANDARDIZE3_PROFILE_TARGET) $(ASSIGN_PERMUTED_SORT_COMPARE_TARGET) $(GRAPH_STANDARDIZER_COMPARE_TARGET) $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty $(RANDOMIZE_SPLIT_CONTRACT_TARGET) $(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(SOLVER_COMPARE_TARGET) $(SOLVER_COMPARE_TARGET)-linbox $(TEST_TARGET)

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

graph-standardizer-compare: $(GRAPH_STANDARDIZER_COMPARE_TARGET)

graph-standardizer-compare-nauty: $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty

run-graph-standardizer-compare: $(GRAPH_STANDARDIZER_COMPARE_TARGET)
	@./$(GRAPH_STANDARDIZER_COMPARE_TARGET) $(or $(WHEEL),25) $(or $(ROUNDS),2) $(or $(REPEAT),1) $(or $(ITER),3) $(or $(WORKLOAD),split-contract)

wheel-split-contract-reps: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)

run-wheel-split-contract-reps: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(or $(WHEEL),11) $(or $(ROUNDS),2) $(OUT)

run-split-contract-files: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "Generating split-contract files for W$(WHEEL), rounds $(SPLIT_CONTRACT_ROUNDS)"
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) $(SPLIT_CONTRACT_CANDIDATES) unused enumerate
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) $(SPLIT_CONTRACT_CANDIDATES) $(SPLIT_CONTRACT_MAP_PREFIX) split-map

run-split-contract-solver: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@echo "Solving from generated split-contract files for W$(WHEEL), rounds $(SPLIT_CONTRACT_ROUNDS)"
	@echo "Checkpoint prefix: $(SPLIT_CONTRACT_MAP_PREFIX)_wiedemann_checkpoint.bin"
	@echo "Checkpoint interval: $(CHECKPOINT_INTERVAL) MMT iterations"
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(WHEEL) $(SPLIT_CONTRACT_ROUNDS) unused $(SPLIT_CONTRACT_MAP_PREFIX) split-map-solve $(CHECKPOINT_INTERVAL)

solver-compare: $(SOLVER_COMPARE_TARGET)

solver-compare-linbox: $(SOLVER_COMPARE_TARGET)-linbox

run-solver-compare: $(SOLVER_COMPARE_TARGET)
	@./$(SOLVER_COMPARE_TARGET) $(or $(ROWS),5000) $(or $(COLS),5000) $(or $(NNZ),50000) $(or $(SEED),17)

test: $(TEST_TARGET)
