# FORAlign: Accelerating gap-affine DNA pairwise sequence alignment using t-blocks based on FOur Russians approach with linear space complexity

FORAlign is a library written in C++17 (C++20 for some features) for speeding up Hirschberg algorithm using Four-Russians approach in multithreads. It runs on Linux and Windows.

## Test Environment 

In order to test our methods, we use the benchmark program developed by [WFA2](https://github.com/smarco/WFA2-lib), and we developed two programs to test our two methods: the `align_benchmark` program has tested `russians-multi` methods, and the `multiswg_benchmark` has tested the `condition-variable` method. The usage for `align_benchmark` is shown as follows, our methods is in \[Gap-affine (Smith-Waterman-Gotoh)\], which called `gap-affine-hirschberg-multi` and `gap-affine-russians-multi`:

```bash
USE: ./align_benchmark -a ALGORITHM -i PATH
      Options::
        [Algorithm]
          --algorithm|a ALGORITHM
            [Indel (Longest Common Subsequence)]
              indel-wfa
            [Edit (Levenshtein)]
              edit-bpm
              edit-dp
              edit-dp-banded
              edit-wfa
            [Gap-linear (Needleman-Wunsch)]
              gap-linear-nw
              gap-linear-wfa
            [Gap-affine (Smith-Waterman-Gotoh)]
              gap-affine-swg
              gap-affine-swg-banded
              gap-affine-wfa
              gap-affine-hirschberg-multi
              gap-affine-russians-multi
            [Gap-affine-2pieces (Concave 2-pieces)]
              gap-affine2p-dp
              gap-affine2p-wfa
        [Input & Output]
          --input|i PATH
          --output|o PATH
          --output-full PATH
          --input-ref|I PATH (Reference alignment file)
        [Penalties]
          --linear-penalties|p M,X,I
          --affine-penalties|g M,X,O,E
          --affine2p-penalties M,X,O1,E1,O2,E2
        [Wavefront parameters]
          --wfa-score-only
          --wfa-span 'global'|'extension'|'ends-free[,P0,Pf,T0,Tf]'
          --wfa-memory 'high'|'med'|'low'|'ultralow'
          --wfa-heuristic STRATEGY
          --wfa-heuristic-parameters  P1,P2[,P3]
            [STRATEGY='banded-static']
              P1 = minimum-diagonal-band (e.g., -100)
              P2 = maximum-diagonal-band (e.g., +100)
            [STRATEGY='banded-adaptive']
              P1 = minimum-diagonal-band (e.g., -100)
              P2 = maximum-diagonal-band (e.g., +100)
              P3 = steps-between-cutoffs
            [STRATEGY='wfa-adaptive']
              P1 = minimum-wavefront-length
              P2 = maximum-difference-distance
              P3 = steps-between-cutoffs
            [STRATEGY='xdrop']
              P1 = x-drop
              P2 = steps-between-cutoffs
            [STRATEGY='zdrop']
              P1 = z-drop
              P2 = steps-between-cutoffs
          --wfa-max-memory BYTES
          --wfa-max-steps INT
          --wfa-max-threads INT (intra-parallelism; default=1)
        [Russians Parameter]
          --russian-block|B INT (Russians block size; default=50)
        [Multithread Parameters]
          --dp-threads|D INT
             (SWG multithread/Hirschberg/Russians DP; default=1)
          --divide-threads|E INT (Hirschberg/Russians main; default=1)
        [Other Parameters]
          --bandwidth INT
        [Misc]
          --check|c 'correct'|'score'|'alignment'
          --check-distance 'indel'|'edit'|'linear'|'affine'|'affine2p'
          --check-bandwidth INT
          --plot
        [System]
          --num-threads|t INT
          --batch-size INT
          --progress|P INT
          --verbose|v INT
          --quiet|q
          --help|h
```

The usage for `multiswg_benchmark` is shown as follows, our method is `gap-affine-swg-multithread`:

```bash
USE: ./align_benchmark -a ALGORITHM -i PATH
      Options::
        [Algorithm]
          --algorithm|a ALGORITHM
            [Gap-affine (Smith-Waterman-Gotoh)]
              gap-affine-swg
              gap-affine-swg-multithread
        [Input & Output]
          --input|i PATH
          --output|o PATH
          --output-full PATH
          --input-ref|I PATH (Reference alignment file)
        [Penalties]
          --affine-penalties|g M,X,O,E
        [Multithread Parameters]
          --dp-threads|D INT
             (SWG multithread DP; default=1)
          --use-barrier|B  use barrier to align (default=false)
        [Misc]
          --check|c 'correct'|'score'|'alignment'
        [System]
          --num-threads|t INT
          --batch-size INT
          --progress|P INT
          --verbose|v INT
          --quiet|q
          --help|h
```

The tests in our paper is shown as follows:

|Test method Name|Test Command<br/>(Run the following programs in `bin` folder)||
|:-:|:-:|:-:|
|benchmark|`./multiswg_benchmark -a gap-affine-swg -i $file_name -o $answer_name`|benchmark will output alignment CIGAR result for testing other methods|
|wfa-high|`./align_benchmark -a gap-affine-wfa --wfa-memory high -i $file_name -I $answer_name -c alignment --wfa-max-threads $cpus`||
|wfa-med|`./align_benchmark -a gap-affine-wfa --wfa-memory med -i $file_name -I $answer_name -c alignment --wfa-max-threads $cpus`||
|wfa-low|`./align_benchmark -a gap-affine-wfa --wfa-memory low -i $file_name -I $answer_name -c alignment --wfa-max-threads $cpus`||
|wfa-ultra-low|`./align_benchmark -a gap-affine-wfa --wfa-memory ultralow -i $file_name -I $answer_name -c alignment --wfa-max-threads $cpus`||
|swg-barrier|`./multiswg_benchmark -a gap-affine-swg-multithread -i $file_name -I $answer_name -B -c alignment -D $cpus`|use `-B` to align with barrier but only supported in `C++20` or newer|
|swg-condition-variable|`./multiswg_benchmark -a gap-affine-swg-multithread -i $file_name -I $answer_name -c alignment -D $cpus`||
|hirschberg-single|`./align_benchmark -a gap-affine-hirschberg-multi -i $file_name -I $answer_name -c alignment -E $cpus -D 1`||
|hirschberg-multi|`./align_benchmark -a gap-affine-hirschberg-multi -i $file_name -I $answer_name -c alignment -E $cpus -D $dpcpus`|`$dpcpus` is same as `expr $cpus / 2`|
|russians-multi-$t$|`./align_benchmark -a gap-affine-russians-multi -i $file_name -I $answer_name -c alignment -E $cpus -D $dpcpus -B $blocksize`|`$dpcpus` is same as `expr $cpus / 2`<br/>`$blocksize` is same as $t$|


## Test dataset and program

The test dataset is stored at [`https://github.com/malabz/FORAlign_testcase`](https://github.com/malabz/FORAlign_testcase). The test cases in paper is stored at [`https://github.com/malabz/FORAlign_testcase/blob/main/test-result/raw.tar.xz`](https://github.com/malabz/FORAlign_testcase/blob/main/test-result/raw.tar.xz), the compiled program is stored at [`https://github.com/malabz/FORAlign_testcase/tree/main/wfa2-test-prog`](https://github.com/malabz/FORAlign_testcase/tree/main/wfa2-test-prog) folder.

### How to run our test program

The WFA2 program is only supported in Linux. For Windows user please install WSL first. WSL instructional video: [1](https://www.youtube.com/watch?v=X-DHaQLrBi8&t=5s) or [2](http://lab.malab.cn/%7Etfr/1.mp4) (Copyright belongs to the original work).

#### Use the compiled WFA2 program

Clone the `FORAlign-testcase` repository and test:
```bash
#1 Download
git clone https://github.com/malabz/FORAlign-testcase

#2 Open the folder
cd wfa2-test-prog

#3 Test
./multiswg_benchmark
./align_benchmark
```

#### Compile WFA2

Clone the `FORAlign` repository:
```bash
#1 Download
git clone https://github.com/malabz/FORAlign

#2 Open the folder and compile
cd FORAlign
make all THREADS=16

#3 Test
./bin/multiswg_benchmark
./bin/align_benchmark

```

## For development

FORAlign work as a programming library. This section shows how to use the C/CPP APIs of FORAlign to take two sequences as input and perform pairwise sequence alignment. Basically, the library file `libforalign.a` and header file `hirschberg.h` are needed to make the Hirschberg algorithm in FORAlign library work in your program, and the library file `libswg.a` and header file `parallel_swg.h` make the SWG algorithm in FORAlign library work in your program.

### A simple example for `libforalign.a`

First, include the FORAlign alignment headers.
```c
#include "foralign/hirschberg.h"
```

Next, configure the function arugments. The hirschberg library uses gap affine penalties. Note that mismatch, gap open penalty, gap extension penalty and russian block size must be positive values, the match penalty should be 0.
```c
/*
 M: match
 X: mismatch
 O: gap open
 E: gap extension
 table_size: russian table size
 threads: program threads
 */
int O, E, M, X, table_size, threads;
threads = 16; O = 6; E = 2; M = 0; X = 4; table_size = 200;
```

Finally, call the `Hirschberg_API` function. If you only want to get CIGAR string, call the `hirschberg_cigar` function.
```c
char *seq1 = "ATCGTAT";
char *seq2 = "TTTTCTAAA";
sequence_t type = DNA;
char *comp_seq1 = NULL, *comp_seq2 = NULL;
hirschberg_API(seq1, seq2, strlen(seq1), strlen(seq2), O, E, M, X, table_size, type, threads, &comp_seq1, &comp_seq2, 0, 1, threads / 2);
```

See `library_example/foralign.c` for C API, `library_example/foralign.cpp` for CPP API. For example, to compile and run example files, you need to link against the FORAlign library (`-lforalign`) as follows:

```bash
make THREADS=16 # make the foralign library
mkdir -p bin
g++ -O3 library_example/foralign.cpp -pthread -L./lib -lforalign -o bin/foralign_cpp_example
# determine the link order
gcc -O3 library_example/foralign.c -o bin/foralign_c_example -L./lib -lforalign -pthread -lstdc++

./bin/foralign_cpp_example
./bin/foralign_c_example
```

### A simple example for `libswg.a`

First, include the FORAlign alignment headers.

```c
#include "swg/parallel_swg.h"
```

Next, configure the function arugments. The swg library uses gap affine penalties. The calcluate matrix should be allocated before call the function:

```c
affine_penalties_t args; args.match = 0; args.mismatch = 4; args.gap_opening = 6; args.gap_extension = 2;
affine_matrix_t    mtx;
mtx.num_columns = len1; mtx.num_rows = len2;
mtx.columns = (affine_cell_t**)malloc(len1 + 1);
if(mtx.columns == nullptr) { fprintf(stderr, "Error: can not allocate. Program will exit.\n"); return 1; }
for(size_t i = 0; i < len2; ++ i)
{
    mtx.columns[i] = (affine_cell_t*)malloc(len2 + 1);
    if(mtx.columns[i] == nullptr) { fprintf(stderr, "Error: can not allocate. Program will exit.\n"); return 1; }
}
```

Finally, call the `multithread_swg_compute_cv` function.

```C
char* seq1 = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
char* seq2 = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
int threads = 16;
multithread_swg_init(threads);
multithread_swg_compute_cv(threads, &mtx, &args, seq1, strlen(seq1), seq2, strlen(seq2));
multithread_swg_free();
```

Afterwards, we can use the library to display the alignment result (e.g., the pairwise alignment result).
```C
char *comp_seq1 = NULL, *comp_seq2 = NULL;
multithread_swg_traceback(&mtx, &args, strlen(seq1), strlen(seq2), seq1, seq2, &comp_seq1, &comp_seq2);
printf("Alignment result: %s\n%s\n", comp_seq1, comp_seq2);
```

See `library_example/swg.cpp` for CPP API. For example, to compile and run example files, you need to link against the FORAlign SWG library (`-lswg`) as follows:

```bash
make THREADS=16 # make the swg and foralign library
mkdir -p bin
g++ -O3 library_example/swg.cpp -pthread -L./lib -lswg -o bin/swg_cpp_example -lforalign
# The function `read_2_seqs` is defined in foralign library

./bin/foralign_swg_cpp_example
```

## Citation

## Contacts

If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/FORAlign/issues) or [email us](mailto:wym6912@outlook.com).

More tools and infomation can visit our [github](https://github.com/malabz).