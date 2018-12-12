#include "seqmini.hpp"
#include <string>


using namespace std;
int main(int argc, char** argv){

    sqm::sqm_pair_list_t* sqt = sqm::minimizers_from_fasta(argv[1]);
    sqm::emit_minimizers_to_text(sqt, "outfile.txt");
    sqm::emit_minimizers_to_binary(sqt, "binary_outfile.txt");
    delete [] sqt;
    return 0;
};