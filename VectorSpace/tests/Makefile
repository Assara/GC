# VectorSpace/tests/Makefile
# Header-only test build: compiles test_*.cpp in this directory into ./bin/
# No src/, no library.

CXX ?= g++
STD ?= -std=c++23

THIS_MAKEFILE := $(lastword $(MAKEFILE_LIST))
THIS_DIR      := $(strip $(abspath $(dir $(THIS_MAKEFILE))))

# VectorSpace/ directory and repo root
VSPACE_ROOT := $(strip $(abspath $(THIS_DIR)/..))
REPO_ROOT   := $(strip $(abspath $(THIS_DIR)/../..))

# Include directories:
#  - repo/include if it exists (common)
#  - repo root (so includes like "VectorSpace/..." work if VectorSpace/ is under repo root)
#  - VectorSpace/ itself (so includes like "lil_matrix.hpp" work if headers live there)
CPPFLAGS += -I"$(REPO_ROOT)/include" -I"$(REPO_ROOT)" -I"$(VSPACE_ROOT)"

CXXFLAGS += $(STD) -O2 -g -Wall -Wextra -Wpedantic -fopenmp
LDFLAGS  += -fopenmp
LDLIBS   +=

BIN_DIR := $(THIS_DIR)/bin

TEST_SRCS := $(wildcard $(THIS_DIR)/test_*.cpp)
TEST_BINS := $(patsubst $(THIS_DIR)/%.cpp,$(BIN_DIR)/%,$(TEST_SRCS))

.PHONY: all run clean rebuild info asan ubsan

all: info $(BIN_DIR) $(TEST_BINS)

info:
	@echo "THIS_DIR     = $(THIS_DIR)"
	@echo "VSPACE_ROOT  = $(VSPACE_ROOT)"
	@echo "REPO_ROOT    = $(REPO_ROOT)"
	@echo "BIN_DIR      = $(BIN_DIR)"
	@echo "TEST_SRCS    = $(notdir $(TEST_SRCS))"
	@echo "TEST_BINS    = $(notdir $(TEST_BINS))"
	@echo "CPPFLAGS     = $(CPPFLAGS)"

$(BIN_DIR):
	@mkdir -p "$(BIN_DIR)"

# Build each test executable directly
$(BIN_DIR)/%: $(THIS_DIR)/%.cpp | $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) "$<" $(LDFLAGS) $(LDLIBS) -o "$@"

run: all
	@set -e; \
	for t in $(TEST_BINS); do \
		echo "==> $$t"; \
		"$$t"; \
	done

asan: CXXFLAGS += -O1 -g -fsanitize=address -fno-omit-frame-pointer
asan: LDFLAGS  += -fsanitize=address
asan: clean all

ubsan: CXXFLAGS += -O1 -g -fsanitize=undefined -fno-omit-frame-pointer
ubsan: LDFLAGS  += -fsanitize=undefined
ubsan: clean all

rebuild: clean all

clean:
	@echo "Cleaning: $(BIN_DIR) and generated binaries"
	rm -rf "$(BIN_DIR)"
