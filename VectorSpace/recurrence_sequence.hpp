#pragma once
#include "packed_recurrence.hpp"
namespace VectorSpace::block_wiedemann_detail {
template<Field K>
void save_recurrence_sequence(const std::filesystem::path& path,const packed_moments<K>& moments,std::size_t training) {
    if(!training || training>moments.size())throw std::invalid_argument("invalid recurrence training length");
    serialization::archive_output out(path);
    out.section("recurrence-comparison",[&]{out.word(1);out.word(moments.block_size());out.word(moments.size());out.word(training);});
    moments.save(out);out.finish();
}
template<Field K>
struct recurrence_sequence {
    packed_moments<K> moments;
    std::size_t training;
    static recurrence_sequence load(const std::filesystem::path& path) {
        serialization::archive_input in(path);
        std::size_t b=0,total=0,training=0;
        in.section("recurrence-comparison",[&]{in.expect(1);b=in.word();total=in.word();training=in.word();});
        if(!b || !training || training>total)throw std::runtime_error("invalid saved recurrence sequence");
        recurrence_sequence result{packed_moments<K>(b,total),training};
        result.moments.load(in);in.finish();return result;
    }
};
}
