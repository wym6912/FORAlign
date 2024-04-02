#include "../foralign/hirschberg.h"
#include <string.h>

int O, E, M, X, table_size;
int threads;
char *seq1 = NULL, *seq2 = NULL, *name1 = NULL, *name2 = NULL;

int main(int argc, char** argv)
{
    threads = 16; O = 6; E = 2; M = 0; X = 4; table_size = 200;
    seq1 = "ATCGTAT";
    seq2 = "TTTTCTAAA";
    char *comp_seq1 = NULL, *comp_seq2 = NULL;
    sequence_t type = DNA;
    // call Hirschberg
    hirschberg_API(seq1, seq2, strlen(seq1), strlen(seq2), O, E, M, X, table_size, type, threads, &comp_seq1, &comp_seq2, 0, 1, threads / 2);
    fprintf(stdout, "Alignment result:\n%s\n%s\n", comp_seq1, comp_seq2);
	return 0;
}
