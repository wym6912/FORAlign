#include "hirschberg.h"
#include "mtxutl.hpp"
#include "../include/kseq/kseq.h"
#include <cstdio>
#include <string>

#if (defined(__linux__) || defined(__APPLE__))
#include <unistd.h>
#else
#include <io.h>
#include <process.h>
#endif

KSEQ_INIT(int, read)

DLL void read_2_seqs(char *c_input, char **seq1, char **name1, char **seq2, char **name2, size_t *len1, size_t *len2)
{
    FILE* f_pointer = fopen(c_input, "r");
    if (f_pointer == NULL) { fprintf(stderr, "Error: file %s cannot open. Program will exit\n", c_input); exit(1); }
    kseq_t* file_t_ = kseq_init(fileno(f_pointer));
    int seq_cnt = 0;
    while (kseq_read(file_t_) >= 0)
    {
        seq_cnt ++;
        if (seq_cnt >= 3)
        {
            fprintf(stderr, "Error: file %s has more than 2 sequences. Program will exit.\n", c_input);
            exit(1);
        }
        if (seq_cnt == 1)
        {
            *seq1 = strdup(file_t_ -> seq.s);
            *len1 = file_t_ -> seq.l;
            if (*seq1 == NULL)
            {
                fprintf(stderr, "Error: Can not save sequence 1. Program will exit.\n");
                exit(1);
            }
            *name1 = AllocateCharVec(file_t_ -> name.l + file_t_ -> comment.l + 1);
            strncpy(*name1, file_t_ -> name.s, file_t_ -> name.l);
            strncpy(*name1 + file_t_ -> name.l, file_t_ -> comment.s, file_t_ -> comment.l);
        }
        else if (seq_cnt == 2)
        {
            *seq2 = strdup(file_t_ -> seq.s);
            *len2 = file_t_ -> seq.l;
            if (*seq2 == NULL)
            {
                fprintf(stderr, "Error: Can not save sequence 2. Program will exit.\n");
                exit(1);
            }
            *name2 = AllocateCharVec(file_t_ -> name.l + file_t_ -> comment.l + 1);
            strncpy(*name2, file_t_ -> name.s, file_t_ -> name.l);
            strncpy(*name2 + file_t_ -> name.l, file_t_ -> comment.s, file_t_ -> comment.l);
        }
    }
    kseq_destroy(file_t_);
    fclose(f_pointer);
}

DLL void print_seqs_to_file(char *seq1, char *comment1, char *seq2, char *comment2, char* c_output)
{
    FILE *f = fopen(c_output, "wb");
    if(f == NULL) { fprintf(stderr, "Error: can not open file %s. Program will exit.\n", c_output); exit(1); }
    fprintf(f, ">%s\n", comment1); fputs(seq1, f); fputc('\n', f);
    fprintf(f, ">%s\n", comment2); fputs(seq2, f); fputc('\n', f);
    fclose(f);
}
