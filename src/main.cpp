#include "seqmini.hpp"
#include <string>


using namespace std;

void usage(char** argv){
    cerr << "seqmini: generate minimizers for arbitrary sequencing files." << endl;
    cerr << "Usage: " << argv[0] << " [options] infile" << endl <<
    "Options: " << endl <<
    "  -f / --fasta   input is in FASTA format." << endl <<
    "  -q / --fastq   input is in FASTQ format." << endl <<
    "  -g / --gfa     input is in GFA format." << endl <<
    "  -s / --seqs-only only create minimizers for GFA S lines (i.e. don't consider graph topology)." << endl <<
    "  -b / --binary  output binary, rather than text." << endl <<
    endl;
};

int main(int argc, char** argv){

    sqm::sqm_pair_list_t* sqt = sqm::minimizers_from_fasta(argv[1]);
    sqm::emit_minimizers_to_text(sqt, "outfile.txt");
    sqm::emit_minimizers_to_binary(sqt, "binary_outfile.txt");
    delete [] sqt;
    return 0;
};