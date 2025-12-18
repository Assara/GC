# ======= CONFIGURATION =======

CXX       := clang++
LD        := clang++

TARGET    := gc
BUILD_DIR := build

# Include directories
INC       := -I. -IVectorSpace

# Compiler and linker flags
CXXFLAGS  := -std=c++23 -O3 -march=native -Wall -Wextra -Wpedantic $(INC) -flto
LDFLAGS   := -flto -fuse-ld=lld

# ======= SOURCE / BUILD SETUP =======

# All .cpp files EXCEPT anything under VectorSpace/tests/
SRCS := $(shell find . -type f -name "*.cpp" -not -path "./VectorSpace/tests/*")

# Map sources to corresponding build objects
OBJS := $(patsubst ./%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# ======= RULES =======

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "🚧 Linking $(TARGET) with Clang + LLD..."
	$(LD) $(OBJS) $(LDFLAGS) -o $(TARGET)
	@echo "✅ Build complete."

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	@echo "🔧 Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@echo "🧹 Cleaning build files..."
	rm -rf $(BUILD_DIR) $(TARGET)

run: all
	@./$(TARGET)
