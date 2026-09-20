#pragma once
#include "mmap_archive.hpp"
#include "../types.hpp"
namespace VectorSpace::block_wiedemann_detail {
// Solver-specific metadata layered over the reusable mmap archive format.
class checkpoint_output : public serialization::archive_output {
public:
    explicit checkpoint_output(const std::filesystem::path& path):archive_output(path) {
        section("solver-types",[&]{
            word(1);word(sizeof(SmallSignedInt));word(sizeof(std::uint32_t));
            word(sizeof(std::size_t));word(sizeof(Int));word(sizeof(GraphAccumulator));
        });
    }
};
class checkpoint_input : public serialization::archive_input {
public:
    explicit checkpoint_input(const std::filesystem::path& path):archive_input(path) {
        section("solver-types",[&]{
            expect(1);expect(sizeof(SmallSignedInt));expect(sizeof(std::uint32_t));
            expect(sizeof(std::size_t));expect(sizeof(Int));expect(sizeof(GraphAccumulator));
        });
    }
};
}
