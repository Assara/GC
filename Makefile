# ======= CONFIGURATION =======

CXX       := clang++
LD        := clang++

TARGET    := gc
PLAYGROUND_TARGET := gc_contraction_playground
SPLIT_WATERFALL_TARGET := gc_split_waterfall_playground
STANDARDIZE_PERF_TARGET := standardize_perf
STANDARDIZE3_PROFILE_TARGET := standardize3_profile
GRAPH_STANDARDIZER_COMPARE_TARGET := graph_standardizer_compare
WHEEL_SPLIT_CONTRACT_REPS_TARGET := wheel_split_contract_reps
SOLVER_COMPARE_TARGET := solver_compare
TEST_TARGET := gc_test
BUILD_DIR := build

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

GRAPH_STANDARDIZER_COMPARE_SRCS := tools/graph_standardizer_comparison.cpp
GRAPH_STANDARDIZER_COMPARE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(GRAPH_STANDARDIZER_COMPARE_SRCS))
GRAPH_STANDARDIZER_COMPARE_NAUTY_OBJS := $(BUILD_DIR)/tools/graph_standardizer_comparison_nauty.o

WHEEL_SPLIT_CONTRACT_REPS_SRCS := tools/wheel_split_contract_reps.cpp
WHEEL_SPLIT_CONTRACT_REPS_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(WHEEL_SPLIT_CONTRACT_REPS_SRCS))

SOLVER_COMPARE_SRCS := tools/solver_comparison.cpp
SOLVER_COMPARE_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOLVER_COMPARE_SRCS))
SOLVER_COMPARE_LINBOX_OBJS := $(BUILD_DIR)/tools/solver_comparison_linbox.o

TEST_SRCS := GC_test.cpp
TEST_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SRCS))

# ======= RULES =======

.PHONY: all clean run playground split-waterfall run-split-waterfall standardize-perf standardize3-profile run-standardize-perf run-standardize3-profile graph-standardizer-compare graph-standardizer-compare-nauty run-graph-standardizer-compare wheel-split-contract-reps run-wheel-split-contract-reps solver-compare solver-compare-linbox run-solver-compare test

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

$(GRAPH_STANDARDIZER_COMPARE_TARGET): $(GRAPH_STANDARDIZER_COMPARE_OBJS)
	@echo "🚧 Linking $(GRAPH_STANDARDIZER_COMPARE_TARGET) with Clang + LLD..."
	$(LD) $(GRAPH_STANDARDIZER_COMPARE_OBJS) $(LDFLAGS) -o $(GRAPH_STANDARDIZER_COMPARE_TARGET)
	@echo "✅ Build complete."

$(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty: $(GRAPH_STANDARDIZER_COMPARE_NAUTY_OBJS)
	@echo "🚧 Linking $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty with Clang + LLD..."
	$(LD) $(GRAPH_STANDARDIZER_COMPARE_NAUTY_OBJS) $(LDFLAGS) $(NAUTY_LDLIBS) -o $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty
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
	rm -rf $(BUILD_DIR) $(TARGET) $(PLAYGROUND_TARGET) $(SPLIT_WATERFALL_TARGET) $(STANDARDIZE_PERF_TARGET) $(GRAPH_STANDARDIZER_COMPARE_TARGET) $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty $(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(SOLVER_COMPARE_TARGET) $(SOLVER_COMPARE_TARGET)-linbox $(TEST_TARGET)

run: all
	@./$(TARGET)

playground: $(PLAYGROUND_TARGET)

split-waterfall: $(SPLIT_WATERFALL_TARGET)

run-split-waterfall: $(SPLIT_WATERFALL_TARGET)
	@./$(SPLIT_WATERFALL_TARGET) $(or $(WHEEL),7) $(or $(OUT),W$(or $(WHEEL),7).txt)

standardize-perf: $(STANDARDIZE_PERF_TARGET)

standardize3-profile: $(STANDARDIZE3_PROFILE_TARGET)

run-standardize-perf: $(STANDARDIZE_PERF_TARGET)
	@./$(STANDARDIZE_PERF_TARGET) $(or $(WHEEL),33) $(or $(REPEAT),1) $(or $(ITER),3)

run-standardize3-profile: $(STANDARDIZE3_PROFILE_TARGET)
	@./$(STANDARDIZE3_PROFILE_TARGET) $(or $(WHEEL),11) $(or $(ROUNDS),2) $(or $(REPEAT),1) $(or $(ITER),3) $(or $(WORKLOAD),split-contract)

graph-standardizer-compare: $(GRAPH_STANDARDIZER_COMPARE_TARGET)

graph-standardizer-compare-nauty: $(GRAPH_STANDARDIZER_COMPARE_TARGET)-nauty

run-graph-standardizer-compare: $(GRAPH_STANDARDIZER_COMPARE_TARGET)
	@./$(GRAPH_STANDARDIZER_COMPARE_TARGET) $(or $(WHEEL),25) $(or $(REPEAT),1) $(or $(ITER),3)

wheel-split-contract-reps: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)

run-wheel-split-contract-reps: $(WHEEL_SPLIT_CONTRACT_REPS_TARGET)
	@./$(WHEEL_SPLIT_CONTRACT_REPS_TARGET) $(or $(WHEEL),11) $(or $(ROUNDS),2) $(OUT)

solver-compare: $(SOLVER_COMPARE_TARGET)

solver-compare-linbox: $(SOLVER_COMPARE_TARGET)-linbox

run-solver-compare: $(SOLVER_COMPARE_TARGET)
	@./$(SOLVER_COMPARE_TARGET) $(or $(ROWS),5000) $(or $(COLS),5000) $(or $(NNZ),50000) $(or $(SEED),17)

test: $(TEST_TARGET)
