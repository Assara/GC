# ======= CONFIGURATION =======

CXX       := clang++
LD        := clang++

TARGET    := gc
PLAYGROUND_TARGET := gc_contraction_playground
TEST_TARGET := gc_test
BUILD_DIR := build

# Include directories
INC       := -I. -IVectorSpace

# ---- OpenMP ----
# For Clang on Linux, libomp is typical. If your system uses libgomp instead, see notes below.
OPENMP_CXXFLAGS := -fopenmp=libomp
OPENMP_LDFLAGS  := -fopenmp=libomp -lomp

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

TEST_SRCS := GC_test.cpp
TEST_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SRCS))

# ======= RULES =======

.PHONY: all clean run playground test

all: $(TARGET)

$(TARGET): $(MAIN_OBJS)
	@echo "🚧 Linking $(TARGET) with Clang + LLD..."
	$(LD) $(MAIN_OBJS) $(LDFLAGS) -o $(TARGET)
	@echo "✅ Build complete."

$(PLAYGROUND_TARGET): $(PLAYGROUND_OBJS)
	@echo "🚧 Linking $(PLAYGROUND_TARGET) with Clang + LLD..."
	$(LD) $(PLAYGROUND_OBJS) $(LDFLAGS) -o $(PLAYGROUND_TARGET)
	@echo "✅ Build complete."

$(TEST_TARGET): $(TEST_OBJS)
	@echo "🚧 Linking $(TEST_TARGET) with Clang + LLD..."
	$(LD) $(TEST_OBJS) $(LDFLAGS) -o $(TEST_TARGET)
	@echo "✅ Build complete."

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo "🧹 Cleaning build files..."
	rm -rf $(BUILD_DIR) $(TARGET) $(PLAYGROUND_TARGET) $(TEST_TARGET)

run: all
	@./$(TARGET)

playground: $(PLAYGROUND_TARGET)

test: $(TEST_TARGET)
