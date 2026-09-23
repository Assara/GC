# Run from the repository root: make -f GraphHomology/Makefile
CXX ?= g++
HOMOLOGY_LOOP ?= 6
HOMOLOGY_VERTICES ?= 8
HOMOLOGY_FLAGS ?= -std=c++23 -O1 -fopenmp -I. -IVectorSpace
BINARY := build/contraction_matrices_L$(HOMOLOGY_LOOP)_V$(HOMOLOGY_VERTICES)
HEADERS := VectorSpace/divide_conquer_generator.hpp VectorSpace/recurrence_sequence.hpp VectorSpace/reconstruction_state.hpp VectorSpace/mmap_archive.hpp VectorSpace/packed_recurrence.hpp VectorSpace/recurrence_checkpoint.hpp GraphHomology/NaturalRepresentatives.hpp GraphHomology/NaturalComposition.hpp GraphHomology/NaturalAdjoint.hpp VectorSpace/OwnedArray.hpp VectorSpace/packed_projection.hpp VectorSpace/sparse_transpose.hpp timer_accum.hpp VectorSpace/sparse_block.hpp GraphHomology/StackedContraction.hpp VectorSpace/block_wiedemann.hpp VectorSpace/block_minimal_generator.hpp VectorSpace/Field.hpp VectorSpace/sparse_rank.hpp VectorSpace/compressed_sparse_matrix.hpp VectorSpace/lil_matrix.hpp VectorSpace/LinComb.hpp VectorSpace/BasisElement.hpp $(wildcard Z*.hpp) GraphHomology/ContractionMatrices.hpp graph.hpp GraphStandardizer.hpp GraphGeneration/CutVertexSplitRule.hpp types.hpp permutation.hpp CombinatorialUtils.hpp graph_hash.hpp
.PHONY: all
all: $(BINARY)
$(BINARY): tools/build_contraction_matrices.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) -DGC_HOMOLOGY_LOOP=$(HOMOLOGY_LOOP) -DGC_HOMOLOGY_VERTICES=$(HOMOLOGY_VERTICES) $< -o $@

build/contraction_matrices_test: tools/contraction_matrices_test.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) -UNDEBUG $< -o $@

.PHONY: test
# Supply the generated L6 graph directory explicitly; no cached output is assumed.
test: build/contraction_matrices_test
	./build/contraction_matrices_test "$(GRAPH_DIRECTORY)"

build/check_natural_adjoint_L$(HOMOLOGY_LOOP)_V$(HOMOLOGY_VERTICES): tools/check_natural_adjoint.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) -DGC_HOMOLOGY_LOOP=$(HOMOLOGY_LOOP) -DGC_HOMOLOGY_VERTICES=$(HOMOLOGY_VERTICES) $< -o $@

build/bench_block_recurrence: tools/bench_block_recurrence.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) $< -o $@

build/packed_recurrence_test: tools/packed_recurrence_test.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) -UNDEBUG $< -o $@

# Standalone recurrence comparison: no graph construction or Krylov timing.
build/bench_recurrence_finders: tools/bench_recurrence_finders.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) $< -o $@

build/recurrence_finders_test: tools/recurrence_finders_test.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) -UNDEBUG $< -o $@

build/find_recurrence: tools/find_recurrence.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) $< -o $@

build/bench_krylov_L$(HOMOLOGY_LOOP)_V$(HOMOLOGY_VERTICES): tools/bench_krylov.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(HOMOLOGY_FLAGS) -DGC_HOMOLOGY_LOOP=$(HOMOLOGY_LOOP) -DGC_HOMOLOGY_VERTICES=$(HOMOLOGY_VERTICES) $< -o $@
