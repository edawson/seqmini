#ifndef SEQMINIHPP
#define SEQMINIHPP
#include <vector>
#include <string>
#include <algorithm>

#include "mkmh.hpp"
#include "gfakluge.hpp"
#include "tinyfa.hpp"


/***
 *  Outline:
 *  SeqMini provides a simple way to index arbitrary
 *  sequences with their minimizers
 *  [ ] FASTA Files
 *  [ ] GFA Files
 *  [ ] vg files
 *  
 *  It produces a single common output format
 *  mapping from minimizer -> list(ID)
 *  This is a binary file
 *  The first 64 bits are the minimizer
 *  The file then alternates 16-bit lengths indicating
 *  the number of bytes in the ID and then the ID in that
 *  number of bytes.
 * 
 *  The text version of this file is similar
 *  the minimizer is written
 *  then a tab
 *  then a comma-separated list of IDs.
 */ 

namespace sqm{

    struct sqm_pair_t{
        char* id;
        mkmh::hash_t min;
        sqm_pair_t(){
            id = NULL;
            min = 0;
        };
        ~sqm_pair_t(){
            //delete [] id;
        };
        sqm_pair_t(char*& id_n, mkmh::hash_t min_n){
            id = id_n;
            min = min_n;
        };
        friend std::ostream& operator<<(std::ostream& os, sqm_pair_t pt){
            os << pt.min << '\t' << pt.id;
            return os;
        };
    };

    bool sqmPairComparator(const sqm_pair_t& a, const sqm_pair_t& b)
    {
        return a.min < b.min;
    };

    struct sqm_pair_list_t{
        uint64_t size;
        uint64_t capacity;
        sqm_pair_t* pairs;

        sqm_pair_list_t(){
            size = 0;
            capacity = 100;
            pairs = new sqm_pair_t[capacity];
        };
        sqm_pair_list_t(int sz){
            size = 0;
            capacity = sz;
            pairs = new sqm_pair_t[capacity];
        };
        ~sqm_pair_list_t(){
            delete [] pairs;
        };
        void resize(double factor = 1.3){
            int newcap = int (this->capacity * factor);
            sqm_pair_t* newpairs = new sqm_pair_t[newcap];
            this->capacity = newcap;
            for (int i = 0; i < this->size; ++i){
                *(newpairs + i) = *(this->pairs + i);
            }
            delete [] this->pairs;
            this->pairs = newpairs;
        };

        void emplace(const sqm_pair_t& p){
            if (size == capacity){
                resize();
            }
            *(this->pairs + size) = p;
            ++size;
        };
    };

    inline void sort_pairs(sqm_pair_list_t*&  pairlist){
        std::sort(pairlist->pairs,
         pairlist->pairs + pairlist->size,
         sqmPairComparator);
    };

    struct sqm_minimizer_index_t{
        uint64_t size = 0;
        uint64_t capacity = 0;
        std::unordered_map<mkmh::hash_t, char*> backing_map;
    };

    inline void emit_minimizers_to_binary(sqm_pair_list_t* pairs, char* outfile){
        ofstream ofi;
        ofi.open(outfile, ios::binary);

        if (ofi.good()){
            for (int i = 0; i < pairs->size; ++i){
                int id_len = strlen(pairs->pairs[i].id);
                ofi.write((char*)&pairs->pairs[i].min, sizeof(mkmh::hash_t));
                ofi.write((char*)&id_len, sizeof(int));
                ofi.write((char*)&pairs->pairs[i].id, id_len * sizeof(char));
            }
        }
        else{
            cerr << "Couldn't open output file " << outfile << "." << endl;
            exit(9);
        }

        ofi.close();
    };
    inline void emit_minimizers_to_text(sqm_pair_list_t* pairs, char* outfile){

        ofstream ofi(outfile);
        if (ofi.good()){
            for (int i = 0; i < pairs->size; ++i){
                ofi << pairs->pairs[i] << endl;
            }
        }
        else{
            cerr << "ERROR: Could not open outfile " << outfile << "." << endl;
        }
        ofi.close();
    };

    inline sqm_pair_list_t* minimizers_from_fasta(const char* filename,
        int k = 16,
        int w = 20,
        bool stream_out = false){

        sqm::sqm_pair_list_t* pair_list = new sqm::sqm_pair_list_t(200);

        TFA::tiny_faidx_t tfi;
        if (!TFA::checkFAIndexFileExists(filename)){
            TFA::createFAIndex(filename, tfi);
            TFA::writeFAIndex(filename, tfi);
        }
        else{
            TFA::parseFAIndex(filename, tfi);
        }

        for (auto s : tfi.seq_to_entry){
            // Since we're iterating an internal map, we don't have to check
            // if the sequence is in said map.
            // Copy the name to a fresh pointer
            char* id = new char[strlen(s.first) + 2];
            id[ s.second->name_len + 1] = '\0';
            strcpy(id, s.first);
            // Copy the sequence to a new pointer
            char* seq;
            TFA::getSequence(tfi, s.first, seq);
            mkmh::mkmh_hash_vec* hv = new mkmh::mkmh_hash_vec(100);
            mkmh::minimizers(seq, s.second->seq_len, k, w, hv);
            //#ifdef DEBUG
            //cerr << hv->size << endl;
            //#endif

            for (int i = 0; i < hv->size; ++i){
                if (stream_out){
                    cout << hv->hashes[i] << '\t' << id << endl;
                }
                else{
                    sqm::sqm_pair_t min_pair(id, hv->hashes[i]);
                    pair_list->emplace(min_pair);
                }

            }

            delete hv;
        }
        
        return pair_list;
    };

    inline sqm_pair_list_t* minimizers_from_gfa(const char* filename, 
    int k = 16, 
    int w = 20,
     bool seqsOnly = true,
      bool emit_to_stream = true){
        gfak::GFAKluge gg;
        gg.parse_gfa_file(filename);
        sqm::sqm_pair_list_t* pair_list = new sqm::sqm_pair_list_t(200);

        // Indexes only the seq lines, without worrying about graph topology
        if (seqsOnly){
            map<string, gfak::sequence_elem, gfak::custom_key> name_to_seq = gg.get_name_to_seq();
            for (auto s : name_to_seq){
                // Since we're iterating an internal map, we don't have to check
                // if the sequence is in said map.
                // Copy the name to a fresh pointer
                char* id = new char[strlen(s.first.c_str()) + 2];
                id[ s.second.name.length() + 1] = '\0';
                strcpy(id, s.first.c_str());
                // Copy the sequence to a new pointer
                char* seq;
                memcpy(seq, s.second.sequence.c_str(), sizeof(char) * s.second.length);
                mkmh::mkmh_hash_vec* hv = new mkmh::mkmh_hash_vec(100);
                mkmh::minimizers(seq, s.second.length, 16, 20, hv);


                for (int i = 0; i < hv->size; ++i){
                    if (emit_to_stream){
                        cout << hv->hashes[i] << '\t' << id << endl;
                    }
                    else{
                        sqm::sqm_pair_t min_pair(id, hv->hashes[i]);
                        pair_list->emplace(min_pair);
                    }

                }
            }
        }
        else{

        }
    return pair_list;
    };

    inline sqm_pair_list_t* minimizers_from_fastq(const char* filename, bool emit_to_stream = true){
        sqm::sqm_pair_list_t* pair_list = new sqm::sqm_pair_list_t(200);

        TFA::tiny_faidx_t tfi;
        if (!TFA::checkFAIndexFileExists(filename)){
            TFA::createFAIndex(filename, tfi);
            TFA::writeFAIndex(filename, tfi);
        }
        else{
            TFA::parseFAIndex(filename, tfi);
        }

        for (auto s : tfi.seq_to_entry){
            // Since we're iterating an internal map, we don't have to check
            // if the sequence is in said map.
            // Copy the name to a fresh pointer
            char* id = new char[strlen(s.first) + 2];
            id[ s.second->name_len + 1] = '\0';
            strcpy(id, s.first);
            // Copy the sequence to a new pointer
            char* seq;
            TFA::getSequence(tfi, s.first, seq);
            mkmh::mkmh_hash_vec* hv = new mkmh::mkmh_hash_vec(100);
            mkmh::minimizers(seq, s.second->seq_len, 21, 30, hv);
            //#ifdef DEBUG
            //cerr << hv->size << endl;
            //#endif

            for (int i = 0; i < hv->size; ++i){
                if (emit_to_stream){
                    cout << hv->hashes[i] << '\t' << id << endl;
                }
                else{
                    sqm::sqm_pair_t min_pair(id, hv->hashes[i]);
                    pair_list->emplace(min_pair);
                }

            }

            delete hv;
        }
        return pair_list;
    };

    inline sqm_pair_list_t* read_binary_pair_file(char* filename);

    inline sqm_pair_list_t* read_text_pair_file(char* filename);




 

}
#endif