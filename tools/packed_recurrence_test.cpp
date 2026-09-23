#include <cassert>
#include <iostream>
#include <fstream>
#include "GraphHomology/ContractionMatrices.hpp"
#include "VectorSpace/block_wiedemann.hpp"
using K=fieldType;
using namespace VectorSpace::block_wiedemann_detail;
void fill(packed_moments<K>& moments,std::size_t end) {
    while(moments.size()<end) {
        auto t=moments.size();auto moment=moments.append();
        for(std::size_t a=0;a<7;++a) {
            K power=1;for(std::size_t i=0;i<t;++i)power*=K(a+1);
            for(std::size_t i=0;i<2;++i)for(std::size_t j=0;j<2;++j)
                moment[i*2+j]+=K((a+1)*(i+1)+3)*power*K((a+2)*(j+1)+1);
        }
    }
}
int main() {
    const auto directory=std::filesystem::path("build")/("checkpoint_test_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory);
    using namespace VectorSpace::serialization;
    assert(crc32c("123456789",9)==0xe3069283U);
    assert(~crc32c_software(~0U,reinterpret_cast<const unsigned char*>("123456789"),9)==0xe3069283U);
    const auto matrix_check=[&]<class Coefficient>() {
        using Matrix=GraphHomology::ContractionMatrixStorage<VectorSpace::OwnedArray,Coefficient>;
        Matrix matrix;matrix.rows=3;
        matrix.offsets=VectorSpace::OwnedArray<std::size_t>(3);
        matrix.offsets[1]=1;matrix.offsets[2]=2;
        matrix.row_indices=VectorSpace::OwnedArray<std::uint32_t>(2);matrix.row_indices[1]=2;
        matrix.coefficients=VectorSpace::OwnedArray<Coefficient>(2);
        matrix.coefficients[0]=-2;matrix.coefficients[1]=3;
        const auto file=directory/("matrix"+std::to_string(sizeof(Coefficient))+".bin");
        {archive_output out(file);matrix.save(out);out.finish();}
        archive_input in(file);auto restored=Matrix::load(in);in.finish();
        assert(restored.rows==3 && restored.columns()==2);
        assert(std::ranges::equal(matrix.offsets,restored.offsets));
        assert(std::ranges::equal(matrix.row_indices,restored.row_indices));
        assert(std::ranges::equal(matrix.coefficients,restored.coefficients));
    };
    matrix_check.operator()<std::int8_t>();matrix_check.operator()<std::int32_t>();
    // Correct checksum and payload, but incompatible coefficient width.
    {
        bool rejected=false;
        try {
            archive_input in(directory/"matrix1.bin");
            using Wide=GraphHomology::ContractionMatrixStorage<VectorSpace::OwnedArray,std::int32_t>;
            auto matrix=Wide::load(in);
        } catch(const std::runtime_error&) {rejected=true;}
        assert(rejected);
    }
    const auto path=directory/"working.bin";
    {
        packed_moments<K> moments(2,10);fill(moments,8);
        minimal_generator_state<K> state(moments,2);state.process_up_to(6);
        checkpoint_output out(path);moments.save(out);state.save(out);out.finish();
    } // Release the small allocation before creating the larger one.
    packed_moments<K> restored(2,24);
    minimal_generator_state<K> state(restored,2);
    {checkpoint_input in(path);restored.load(in);state.load(in);in.finish();}
    assert(restored.size()==8 && state.processed_terms()==6);
    const auto bytes=state.allocated_bytes();const auto* data=restored.data();
    fill(restored,24);state.process_up_to(16);
    assert(bytes==state.allocated_bytes() && restored.data()==data);
    assert(restored[1].data()==restored[0].data()+4);
    packed_generator<K> result(2,12);
    assert(state.generator(result));
    assert(result.coefficient(1,0).data()==result.data()+26);
    packed_moments<K> reference(2,24);fill(reference,24);
    auto expected=minimal_generator(reference,2,16);
    assert(expected);
    for(std::size_t r=0;r<2;++r)assert(std::ranges::equal(result.row(r),expected->row(r)));
    {checkpoint_output out(directory/"final.bin");result.save(out);out.finish();}
    packed_generator<K> loaded(2,30);
    {checkpoint_input in(directory/"final.bin");loaded.load(in);in.finish();}
    for(std::size_t r=0;r<2;++r)assert(std::ranges::equal(result.row(r),loaded.row(r)));
    // Damage header metadata, then the checksum, and reject both.
    for(bool metadata:{true,false}) {
        const auto bad=directory/(metadata?"bad_type.bin":"bad_checksum.bin");
        std::filesystem::copy_file(path,bad);
        {std::fstream file(bad,std::ios::binary|std::ios::in|std::ios::out);
            if(metadata)file.seekp(24);else file.seekp(-1,std::ios::end);
            char x=char(0xff);file.write(&x,1);
            if(metadata) {
                file.flush();std::vector<unsigned char> bytes(std::filesystem::file_size(bad));
                file.seekg(0);file.read(reinterpret_cast<char*>(bytes.data()),bytes.size());
                encode_word(bytes.data()+bytes.size()-8,crc32c(bytes.data(),bytes.size()-8));
                file.seekp(0);file.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());
            }}
        bool rejected=false;
        try {checkpoint_input in(bad);restored.load(in);state.load(in);in.finish();}
        catch(const std::runtime_error&){rejected=true;}
        assert(rejected);
    }
    // Valid checksum but wrong field specifically in a per-vector descriptor.
    {
        checkpoint_output out(directory/"bad_vector_field.bin");
        out.section("current-krylov",[&]{out.text(K::name());out.word(K::characteristic());out.word(sizeof(K));
        out.word(1);out.word(1);out.word(1);
        out.word(0);out.word(3);out.word(sizeof(K));out.word(1);
        K value=1;out.bytes(&value,sizeof(value));});out.finish();
    }
    bool rejected=false;
    try {
        checkpoint_input in(directory/"bad_vector_field.bin");K value;
        in.block<K>("current-krylov",std::span<K>(&value,1),1,1);in.finish();
    }catch(const std::runtime_error&){rejected=true;}
    assert(rejected);
    // Exercise the real solver: exhaust capacity, release, resume with more RAM.
    using Solver=VectorSpace::block_wiedemann_solver<K>;
    std::size_t uninterrupted_calls=0,resumed_calls=0;
    auto apply=[&](std::size_t& calls){return [&calls](auto in,auto out,std::size_t b){
        ++calls;for(std::size_t i=0;i<64;++i)for(std::size_t j=0;j<b;++j)
            out[i*b+j]=i<63 ? K(i+1)*in[i*b+j] : K{};
    };};
    for(bool online:{false,true}) {
        Solver::options config{2,1,8,17,1};config.sequence_capacity=140;
        config.incremental_recurrence=online;
        auto full=Solver::from_square_operator(64,apply(uninterrupted_calls),config).nullspace();
        config.sequence_capacity=12;config.checkpoint_path=directory/(online?"online.bin":"batch.bin");
        bool exhausted=false;
        try {auto partial=Solver::from_square_operator(64,apply(resumed_calls),config).nullspace();}
        catch(const std::length_error&){exhausted=true;}
        assert(exhausted && std::filesystem::exists(config.checkpoint_path));
        config.sequence_capacity=140;config.resume_checkpoint=true;
        auto resumed=Solver::from_square_operator(64,apply(resumed_calls),config).nullspace();
        assert(full.complete && resumed.complete && full.basis==resumed.basis);
    }
    assert(resumed_calls==uninterrupted_calls+2); // Only restart's initial B*V is repeated.
    // Interrupt Horner evaluation, including a later attempt with accepted vectors.
    for(std::size_t stop_attempt:{1,2}) {
        Solver::options config{2,1,8,17,1};config.sequence_capacity=140;
        std::size_t complete_calls=0,split_calls=0,horner_calls=0;
        std::ostringstream log;config.log=&log;
        bool interrupt=false;
        const auto snapshot=directory/("reconstruction"+std::to_string(stop_attempt));
        auto operation=[&](std::size_t& calls) {return [&calls,&interrupt,&horner_calls,&log,&snapshot,stop_attempt](auto in,auto out,std::size_t b) {
            if(interrupt && std::filesystem::exists(snapshot.string()+".reconstruction")
                && log.str().find("nullspace attempt "+std::to_string(stop_attempt)+"/")!=std::string::npos
                && log.str().find("Horner reconstruction",log.str().find("nullspace attempt "+std::to_string(stop_attempt)+"/"))!=std::string::npos
                && ++horner_calls==4)throw std::runtime_error("simulated interruption");
            ++calls;
            for(std::size_t i=0;i<64;++i)for(std::size_t j=0;j<b;++j)
                out[i*b+j]=i<59 ? K(i+1)*in[i*b+j] : K{};
        };};
        auto full=Solver::from_square_operator(64,operation(complete_calls),config).nullspace();
        assert(full.complete && full.basis.size()==5);
        log.str("");log.clear();config.checkpoint_path=snapshot;
        config.reconstruction_checkpoint_seconds=0;interrupt=true;
        bool interrupted=false;
        try {auto partial=Solver::from_square_operator(64,operation(split_calls),config).nullspace();}
        catch(const std::runtime_error& error) {
            interrupted=std::string(error.what())=="simulated interruption";
            if(!interrupted)throw;
        }
        assert(interrupted);interrupt=false;config.resume_checkpoint=true;
        const auto prior_log=log.str().size();
        auto resumed=Solver::from_square_operator(64,operation(split_calls),config).nullspace();
        assert(resumed.complete && resumed.basis==full.basis);
        assert(log.str().substr(prior_log).find("reconstruction resumed")!=std::string::npos);
        assert(split_calls==complete_calls+1); // Only the saved-operator check repeats.
        // A completed Horner checkpoint also reloads without an extra Horner step.
        auto again=Solver::from_square_operator(64,operation(split_calls),config).nullspace();
        assert(again.complete && again.basis==full.basis);
    }
    // A rank-only process hands its validated recurrence to a fresh solver.
    for(bool wide:{false,true}) {
        Solver::options config{4,1,8,17,1};config.sequence_capacity=140;
        std::ostringstream log;config.log=&log;
        const auto forward=[](auto in,auto out,std::size_t b) {
            for(std::size_t i=0;i<out.size()/b;++i)for(std::size_t j=0;j<b;++j)
                out[i*b+j]=i<29 ? K(i+1)*in[i*b+j] : K{};
        };
        auto make=[&]{return wide ? Solver(31,32,forward,forward,config)
            : Solver::from_square_operator(32,forward,config);};
        auto full=make().nullspace();assert(full.complete);
        config.checkpoint_path=directory/(wide?"wide-rank":"square-rank");
        auto rank=make().rank();assert(rank.rank==29 && rank.nullity==3);
        assert(std::filesystem::exists(config.checkpoint_path.string()+".rank"));
        config.resume_checkpoint=true;log.str("");log.clear();
        auto solver=make();auto resumed=solver.nullspace();
        assert(resumed.complete && resumed.basis==full.basis);
        assert(log.str().find("rank handoff loaded")!=std::string::npos);
        assert(solver.sequence_stats.sequences==0); // No repeated Krylov or recurrence.
        assert(log.str().find("rank trial")==std::string::npos);
        const auto again=make().rank();assert(again.rank==rank.rank);
        config.seed++;
        bool rejected=false;
        try {auto incompatible=make().rank();}catch(const std::runtime_error&){rejected=true;}
        assert(rejected);
    }
    // Interrupt after durable saves throughout the production sequence path.
    for(const std::string phase:{"krylov","recurrence","validation","validation-complete"}) {
        Solver::options config{2,1,8,17,1};config.sequence_capacity=140;
        std::size_t full_calls=0,split_calls=0,saves=0;
        auto full=Solver::from_square_operator(64,apply(full_calls),config).nullspace();
        config.checkpoint_path=directory/("periodic-"+phase);config.checkpoint_seconds=0;
        config.on_checkpoint=[&](std::string_view saved_phase) {
            if(saved_phase==phase && ++saves==(phase=="validation-complete"?1:3))throw std::runtime_error("checkpoint interruption");
        };
        bool interrupted=false;
        try {auto partial=Solver::from_square_operator(64,apply(split_calls),config).nullspace();}
        catch(const std::runtime_error& e){interrupted=std::string(e.what())=="checkpoint interruption";if(!interrupted)throw;}
        assert(interrupted);config.on_checkpoint={};config.resume_checkpoint=true;
        auto resumed=Solver::from_square_operator(64,apply(split_calls),config).nullspace();
        assert(resumed.complete && resumed.basis==full.basis);
        assert(split_calls==full_calls+1); // Only the starting-block consistency check.
    }
    // A later nullspace block has its own sequence checkpoint and can resume
    // even when an earlier completed reconstruction checkpoint also exists.
    {
        Solver::options config{2,1,8,17,1};config.sequence_capacity=140;
        const auto operation=[](auto in,auto out,std::size_t b){
            for(std::size_t i=0;i<64;++i)for(std::size_t j=0;j<b;++j)
                out[i*b+j]=i<59 ? K(i+1)*in[i*b+j] : K{};
        };
        auto full=Solver::from_square_operator(64,operation,config).nullspace();
        config.checkpoint_path=directory/"later-sequence";config.checkpoint_seconds=0;
        std::ostringstream log;config.log=&log;std::size_t saves=0;
        config.on_checkpoint=[&](std::string_view phase){
            if(phase=="krylov" && log.str().find("nullspace attempt 2/")!=std::string::npos && ++saves==3)
                throw std::runtime_error("later sequence interruption");
        };
        bool interrupted=false;
        try {auto partial=Solver::from_square_operator(64,operation,config).nullspace();}
        catch(const std::runtime_error& e){interrupted=std::string(e.what())=="later sequence interruption";if(!interrupted)throw;}
        assert(interrupted);config.on_checkpoint={};config.resume_checkpoint=true;
        log.str("");log.clear();
        auto resumed=Solver::from_square_operator(64,operation,config).nullspace();
        assert(resumed.complete && resumed.basis==full.basis);
        assert(log.str().find("checkpoint resumed moments=3")!=std::string::npos);
    }
    // Restore in the middle of validation, preserving already checked batches.
    {
        packed_moments<K> moments(2,520);fill(moments,520);
        minimal_generator_state<K> state(moments,2);state.process_up_to(512);
        packed_generator<K> generator(2,256);
        const auto path=directory/"validation-progress.bin";
        bool interrupted=false;
        try {
            state.generator(generator,[&](auto done,auto total){
                if(done && done<total) {
                    checkpoint_output out(path);state.save(out);out.finish();
                    throw std::runtime_error("validation interruption");
                }
            });
        }catch(const std::runtime_error&){interrupted=true;}
        assert(interrupted);
        minimal_generator_state<K> restored(moments,2);
        {checkpoint_input in(path);restored.load(in);in.finish();}
        bool first=true;
        assert(restored.generator(generator,[&](auto done,auto){
            if(done && first){assert(done>256);first=false;}
        }));
        assert(!first);
    }
    std::cout << "Packed storage, larger-capacity restore, type/field/checksum rejection and all-stage periodic resume passed\n";
}
