# ======= CONFIGURATION =======

# Use Clang instead of GCC
CXX       := clang++
LD        := clang++

# SuiteSparse include dir (on Ubuntu)
SUITESPARSE_INC := -I/usr/include/suitesparse

# Explicit include directories (adjust as needed)
INC       := -I. -IVectorSpace -Ioutput $(SUITESPARSE_INC)

# SuiteSparse libs
SUITESPARSE_LIBS := \
    -lspqr -lcholmod -lsuitesparseconfig \
    -lamd -lcolamd -lcamd -lccolamd -ldl

# Compiler and linker flags
CXXFLAGS  := -std=c++23 -O3 -march=native -Wall -Wextra -Wpedantic $(INC) -flto
LDFLAGS   := -flto -fuse-ld=lld $(SUITESPARSE_LIBS)

# Output binary
TARGET    := gc

# ======= SOURCE / BUILD SETUP =======

# Find all .cpp files recursively
SRCS := $(shell find . -type f -name "*.cpp")

# Where to place .o files
BUILD_DIR := build

# Map sources to corresponding build objects
OBJS := $(patsubst ./%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# ======= RULES =======

all: $(TARGET)

# Link everything into one optimized binary
$(TARGET): $(OBJS)
	@echo "🚧 Linking $(TARGET) with Clang + LLD + SuiteSparse..."
	$(LD) $(OBJS) $(LDFLAGS) -o $(TARGET)
	@echo "✅ Build complete."

# Compile each .cpp → .o (mirroring folder structure)
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Remove build artifacts
clean:
	@echo "🧹 Cleaning build files..."
	rm -rf $(BUILD_DIR) $(TARGET)

# Run after build
run: all
	@./$(TARGET)

.PHONY: all clean run
