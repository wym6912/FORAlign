/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "multiswg_benchmark_params.h"

/*
 * Default parameters
 */
align_bench_params_t parameters = {
  // Algorithm
  .algorithm = alignment_gap_affine_swg,
  // I/O
  .input_filename = NULL,
  .output_filename = NULL,
  .input_ref_filename = NULL,
  .output_full = false,
  .output_file = NULL,
  // I/O internals
  .input_file = NULL,
  .line1 = NULL,
  .line2 = NULL,
  .line1_allocated = 0,
  .line2_allocated = 0,
  // Penalties
  .affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  },
  // swg multithread
  .dp_threads = 1,
  .use_barrier = false,
  // Misc
  .check_bandwidth = -1,
  .check_display = false,
  .check_correct = false,
  .check_score = false,
  .check_alignments = false,
  .check_metric = ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE,
  .plot = 0,
  // System
  .num_threads = 1,
  .batch_size = 10000,
  .progress = 100000,
  .verbose = 0,
};
/*
 * Menu
 */
void usage() {
  fprintf(stderr,
      "USE: ./align_benchmark -a ALGORITHM -i PATH                             \n"
      "      Options::                                                         \n"
      "        [Algorithm]                                                     \n"
      "          --algorithm|a ALGORITHM                                       \n"
      "            [Gap-affine (Smith-Waterman-Gotoh)]                         \n"
      "              gap-affine-swg                                            \n"
      "              gap-affine-swg-multithread                                \n"
      "        [Input & Output]                                                \n"
      "          --input|i PATH                                                \n"
      "          --output|o PATH                                               \n"
      "          --output-full PATH                                            \n"
      "          --input-ref|I PATH (Reference alignment file)                 \n"
      "        [Penalties]                                                     \n"
      "          --affine-penalties|g M,X,O,E                                  \n"
      "        [Multithread Parameters]                                        \n"
      "          --dp-threads|D INT                                            \n"
      "             (SWG multithread DP; default=1)                            \n"
      "          --use-barrier|B  use barrier to align (default=false)         \n"
      "        [Misc]                                                          \n"
      "          --check|c 'correct'|'score'|'alignment'                       \n"
      "        [System]                                                        \n"
      "          --num-threads|t INT                                           \n"
      "          --batch-size INT                                              \n"
      "          --progress|P INT                                              \n"
      "          --verbose|v INT                                               \n"
      "          --quiet|q                                                     \n"
      "          --help|h                                                      \n");
}
/*
 * Parsing arguments
 */
void parse_arguments(
    int argc,
    char** argv) {
  struct option long_options[] = {
    /* Input */
    { "algorithm", required_argument, 0, 'a' },
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    { "output-full", required_argument, 0, 800 },
    { "input-ref", required_argument, 0, 'I' },
    /* Penalties */
    { "affine-penalties", required_argument, 0, 'g' },
    /* Multithread swg parameters */
    { "dp-threads", required_argument, 0, 'D' },
    { "use-barrier", no_argument, 0, 'B' },
    /* Misc */
    { "check", required_argument, 0, 'c' },
    /* System */
    { "num-threads", required_argument, 0, 't' },
    { "batch-size", required_argument, 0, 4000 },
    { "progress", required_argument, 0, 'P' },
    { "verbose", required_argument, 0, 4001 },
    { "verbose1", no_argument, 0, 'v' },
    { "quiet", no_argument, 0, 'q' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"a:i:o:p:g:P:c:vqt:hD:E:BI:",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * Algorithm
     */
    case 'a': {
      /* Test bench */
      // Gap-Affine
      if (strcmp(optarg,"gap-affine-swg")==0 ||
          strcmp(optarg,"gap-affine-dp")==0) {
        parameters.algorithm = alignment_gap_affine_swg;
      } else if (strcmp(optarg, "gap-affine-swg-multithread") == 0) {
        parameters.algorithm = alignment_gap_affine_swg_multi;
      } else {
        fprintf(stderr,"Algorithm '%s' not recognized\n",optarg);
        exit(1);
      }
      break;
    }
    /*
     * Input/Output
     */
    case 'i':
      parameters.input_filename = optarg;
      break;
    case 'o':
      parameters.output_filename = optarg;
      break;
    case 800: // --output-full
      parameters.output_filename = optarg;
      parameters.output_full = true;
      break;
    case 'I':
      parameters.input_ref_filename = optarg;
      break;
    /*
     * Penalties
     */
    case 'g': { // --affine-penalties M,X,O,E
      char* sentinel = strtok(optarg,",");
      parameters.affine_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_opening = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_extension = atoi(sentinel);
      break;
    }
    /* 
     * Russians, Hirschberg and Multithread swg parameters 
     */
    case 'D':
      parameters.dp_threads = atoi(optarg);
      break;
    case 'B':
      parameters.use_barrier = true;
      break;
    /*
     * Misc
     */
    case 'c':
      if (optarg ==  NULL) { // default = correct
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"display")==0) {
        parameters.check_display = true;
      } else if (strcasecmp(optarg,"correct")==0) {
        parameters.check_correct = true;
        parameters.check_score = false;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"score")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = false;
      } else if (strcasecmp(optarg,"alignment")==0) {
        parameters.check_correct = true;
        parameters.check_score = true;
        parameters.check_alignments = true;
      } else {
        fprintf(stderr,"Option '--check' must be in {'correct','score','alignment'}\n");
        exit(1);
      }
      break;
    /*
     * System
     */
    case 't': // --num-threads
      parameters.num_threads = atoi(optarg);
      break;
    case 4000: // --batch-size
      parameters.batch_size = atoi(optarg);
      break;
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'v':
      parameters.verbose = 1;
      break;
    case 4001: // --verbose (long option)
      parameters.verbose = atoi(optarg);
      if (parameters.verbose < 0 || parameters.verbose > 4) {
        fprintf(stderr,"Option '--verbose' must be in {0,1,2,3,4}\n");
        exit(1);
      }
      break;
    case 'q':
      parameters.verbose = -1;
      break;
    case 'h':
      usage();
      exit(1);
    // Other
    case '?': default:
      fprintf(stderr,"Option not recognized \n");
      exit(1);
    }
  }
  // Checks general
  if (parameters.input_filename==NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
}
