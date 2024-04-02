#include "../foralign/hirschberg.h"
#include <cstring>

#if (defined(__linux__) || defined(__APPLE__))
#include <getopt.h>
#else
#include "../include/getopt9/include/getopt.h"
#include <string>
#endif

char* c_input = NULL, * c_output = NULL;
int O, E, M, X, table_size;
int threads;
char *seq1 = NULL, *seq2 = NULL, *name1 = NULL, *name2 = NULL;
size_t len1 = 0, len2 = 0;

void version()
{
    fprintf(stderr, "Version 1.0.0.0\n");
    exit(0);
}

void usage(char name[])
{
    fprintf(stderr, "\n\tUsage: %s [options]\n\n", name);
    fprintf(stderr, "Pairwise sequence alignment using Four Russians algorithm speedups parallel Hirschberg with affine gap penalty\n");
    fprintf(stderr, "DP bits = %zu\n", sizeof(dp_t));
    fprintf(stderr, "Available options:\n");
    fprintf(stderr, "\t--in        FILE      sequence file name (Required)\n");
    fprintf(stderr, "\t--out       FILE      output file name (Required)\n");
    fprintf(stderr, "\t--threads   N         use N threads (N >= 1, default: 1), DP will use N/2 threads\n");
    fprintf(stderr, "\t--open      O         gap open penalty (default: 3)\n");
    fprintf(stderr, "\t--ext       E         gap extension penalty (default: 1)\n");
    fprintf(stderr, "\t--match     M         match score (default: 0)\n");
    fprintf(stderr, "\t--mismatch  X         mismatch score (default: 2)\n");
    fprintf(stderr, "\t--block     B         russians block size (default: 200)\n");
    fprintf(stderr, "\t--help                print help message\n");
    fprintf(stderr, "\t--version             show program version\n");
    fprintf(stderr, "Example:\n\t%s --in seq.fasta --out seq_out.fasta\n", name);
}


int main(int argc, char** argv)
{
    threads = 1; O = 6; E = 2; M = 0; X = 4; table_size = 200;
#if ( defined(_WIN32) || defined(_WIN64) )
    std::string consts[] = { "in", "out", "help", "version", "threads", "open", "ext", "match", "mismatch", "block" };
#endif
    int c;

    while (1)
    {
        int option_index = 0;
#if ( defined(_WIN32) || defined(_WIN64) )
        static struct option long_options[] =
        {
            {(char*)consts[0].c_str(), required_argument, 0, 'i'},
            {(char*)consts[1].c_str(), required_argument, 0, 'o'},
            {(char*)consts[2].c_str(), no_argument,       0, 'h'},
            {(char*)consts[3].c_str(), no_argument,       0, 'v'},
            {(char*)consts[4].c_str(), required_argument, 0, 't'},
            {(char*)consts[5].c_str(), required_argument, 0, 'O'},
            {(char*)consts[6].c_str(), required_argument, 0, 'E'},
            {(char*)consts[7].c_str(), required_argument, 0, 'M'},
            {(char*)consts[8].c_str(), required_argument, 0, 'X'},
            {(char*)consts[9].c_str(), required_argument, 0, 'B'},
            {0,                        0,                 0,  0 }
        };
#else
        static struct option long_options[] =
        {
            {"in",       required_argument, 0, 'i'},
            {"out",      required_argument, 0, 'o'},
            {"help",     no_argument,       0, 'h'},
            {"version",  no_argument,       0, 'v'},
            {"threads",  required_argument, 0, 't'},
            {"open",     required_argument, 0, 'O'},
            {"ext",      required_argument, 0, 'E'},
            {"match",    required_argument, 0, 'M'},
            {"mismatch", required_argument, 0, 'X'},
            {"block",    required_argument, 0, 'B'},
            {0,          0,                 0,  0 }
        };
#endif
        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1) break;
        switch (c)
        {
        case 0:
            fprintf(stderr, "Not supported: option %s", long_options[option_index].name);
            if (optarg) fprintf(stderr, " with arg %s", optarg);
            fprintf(stderr, "\n");
            break;
        case 'h':
            usage(argv[0]);
            version();
            break;
        case 'i':
            c_input = argv[optind - 1];
            break;
        case 'o':
            c_output = argv[optind - 1];
            break;
        case 't':
            threads = atoi(argv[optind - 1]);
            break;
        case 'O':
            O = atoi(argv[optind - 1]);
            break;
        case 'E':
            E = atoi(argv[optind - 1]);
            break;
        case 'M':
            M = atoi(argv[optind - 1]);
            break;
        case 'X':
            X = atoi(argv[optind - 1]);
            break;
        case 'B':
            table_size = atoi(argv[optind - 1]);
            break;
        case 'v':
            version();
        }
    }
    if(c_input == NULL || c_output == NULL)
    {
        fprintf(stderr, "Error: not specified input or output file. Check your arugments and try again.\n");
        exit(1);
    }
    // read 2 sequences
    read_2_seqs(c_input, &seq1, &name1, &seq2, &name2, &len1, &len2);
    fprintf(stderr, "Hirschberg parallel\n");
    char *comp_seq1 = nullptr, *comp_seq2 = nullptr;
    sequence_t type = DNA;
    // call Hirschberg
    hirschberg_API(seq1, seq2, strlen(seq1), strlen(seq2), O, E, M, X, table_size, type, threads, &comp_seq1, &comp_seq2, false, true, threads / 2);
    // print sequences
    print_seqs_to_file(comp_seq1, name1, comp_seq2, name2, c_output);
	return 0;
}
