#include "seqmini.hpp"
#include <string>
#include <getopt.h>


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
    "  -w / --window  window size." << endl <<
    "  -k / --kmer    kmer size." << endl <<
    endl;
};

int main(int argc, char** argv){

    string infile;
    bool binary_out = false;
    bool seqs_only = true;
    string input_type = "FASTA";
    int k = 16;
    int w = 20;

    if (argc < 2){
        cerr << "No input given. Please provide a GFA, FASTQ, or FASTA file." << endl;
        usage(argv);
        exit(1);
    }

    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"fasta", no_argument, 0, 'f'},
            {"fastq", no_argument, 0, 'q'},
            {"gfa", no_argument, 0, 'g'},
            {"seqs-only", no_argument, 0, 's'},
            {"kmer", no_argument, 0, 'k'},
            {"window", no_argument, 0, 'w'},
            {"binary", no_argument, 0, 'b'},
            {"version", no_argument, 0, 'v'},
            {0,0,0,0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "k:w:hfqgsbv", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case '?':
            case 'h':
                usage(argv);
                exit(0);
            case 'f':
                input_type = "FASTA";
                break;
            case 'q':
                input_type = "FASTQ";
                break;
            case 'g':
                input_type = "GFA";
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'w':
                w = atoi(optarg);
                break;
            case 's':
                seqs_only = true;
                break;
            case 'b':
                binary_out = true;
                break;
            default:
                abort();
        }
    }

    if (optind > argc){
        cerr << "Error: no input file provided." << endl;
        exit(1);
    }

    infile = argv[optind];

    if (input_type == "FASTA"){
        //sqm::sqm_pair_list_t* sqt = 
        sqm::minimizers_from_fasta(infile.c_str(), k, w, true);
        //sqm::emit_minimizers_to_text(sqt, "outfile.txt");
        //sqm::emit_minimizers_to_binary(sqt, "binary_outfile.txt");
        //delete sqt;
    }
    else if (input_type == "FASTQ"){
        //sqm::sqm_pair_list_t* sqt = 
        sqm::minimizers_from_fasta(infile.c_str(), k, w, true);
        //sqm::emit_minimizers_to_text(sqt, "outfile.txt");
        //sqm::emit_minimizers_to_binary(sqt, "binary_outfile.txt");
        //delete sqt;
    }
    else if (input_type == "GFA"){
        sqm::minimizers_from_gfa(infile.c_str(), k, w, true, true);
    }



    return 0;
};