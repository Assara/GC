#pragma once
#include <limits>
#include <span>
#include <stdexcept>
#include "Field.hpp"
#include "OwnedArray.hpp"
#include "recurrence_checkpoint.hpp"

namespace VectorSpace::block_wiedemann_detail {
inline std::size_t recurrence_product(std::size_t a,std::size_t b) {
    if(b && a>std::numeric_limits<std::size_t>::max()/b)
        throw std::length_error("recurrence allocation size overflow");
    return a*b;
}

// Fixed capacity, allocated at setup. Appending returns the next contiguous block.
template<Field K>
class packed_moments {
    std::size_t block_,used_=0;
    OwnedArray<K> values_;
public:
    packed_moments(std::size_t block,std::size_t capacity)
        :block_(block),values_(recurrence_product(recurrence_product(block,block),capacity)) {
        if(!block)throw std::invalid_argument("zero recurrence block width");
    }
    std::size_t block_size() const {return block_;}
    std::size_t size() const {return used_;}
    std::size_t capacity() const {return values_.size()/(block_*block_);}
    const K* data() const {return values_.data();}
    std::span<K> operator[](std::size_t t) {return {values_.data()+t*block_*block_,block_*block_};}
    std::span<const K> operator[](std::size_t t) const {return {values_.data()+t*block_*block_,block_*block_};}
    std::span<K> append() {
        if(used_==capacity())throw std::length_error("projected sequence capacity exhausted; increase sequence_capacity at setup");
        return (*this)[used_++];
    }
    void save(serialization::archive_output& out) const {
        out.section("projected-coefficients",[&] {
            out.word(block_);out.word(used_);
            out.block<K>("projected-sequence",std::span<const K>(values_.data(),used_*block_*block_),used_*block_,block_);
        });
    }
    void load(serialization::archive_input& in) {
        in.section("projected-coefficients",[&] {
            in.expect(block_);const auto used=in.word();
            if(used>capacity())throw std::length_error("checkpoint needs a larger sequence capacity");
            in.block<K>("projected-sequence",std::span<K>(values_.data(),used*block_*block_),used*block_,block_);
            used_=used;
        });
    }
    std::size_t allocated_bytes() const {return values_.size()*sizeof(K);}
};

// One coefficient slab, one small degree array. Each row has a fixed stride;
// coefficient access is base + row*stride + power*block, with no row pointers.
template<Field K>
class packed_generator {
    OwnedArray<std::size_t> degrees_;
    OwnedArray<K> values_;
public:
    packed_generator(std::size_t block,std::size_t max_degree)
        :degrees_(block),values_(recurrence_product(recurrence_product(block,block),max_degree+1)) {
        if(!block || max_degree==std::numeric_limits<std::size_t>::max())
            throw std::invalid_argument("invalid recurrence dimensions");
    }
    std::size_t size() const {return degrees_.size();}
    std::size_t degree(std::size_t row) const {return degrees_[row];}
    std::size_t max_degree() const {return values_.size()/size()/size()-1;}
    void set_degree(std::size_t row,std::size_t degree) {
        if(degree>max_degree())throw std::length_error("recurrence degree capacity exhausted");
        degrees_[row]=degree;
    }
    std::span<K> coefficient(std::size_t row,std::size_t power) {
        return {values_.data()+row*(max_degree()+1)*size()+power*size(),size()};
    }
    std::span<const K> coefficient(std::size_t row,std::size_t power) const {
        return {values_.data()+row*(max_degree()+1)*size()+power*size(),size()};
    }
    std::span<const K> row(std::size_t r) const {
        return {coefficient(r,0).data(),(degree(r)+1)*size()};
    }
    void save(serialization::archive_output& out) const {
        out.section("final-recurrence",[&] {
            out.word(size());
            for(std::size_t r=0;r<size();++r) {
                out.word(degree(r));out.block<K>("final-recurrence-row",row(r),degree(r)+1,size());
            }
        });
    }
    void load(serialization::archive_input& in) {
        in.section("final-recurrence",[&] {
            in.expect(size());
            for(std::size_t r=0;r<size();++r) {
                const auto d=in.word();set_degree(r,d);
                in.block<K>("final-recurrence-row",std::span<K>(coefficient(r,0).data(),(d+1)*size()),d+1,size());
            }
        });
    }
    const K* data() const {return values_.data();}
    std::size_t allocated_bytes() const {return values_.size()*sizeof(K)+degrees_.size()*sizeof(std::size_t);}
};
}
