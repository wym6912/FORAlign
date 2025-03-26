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
 * DESCRIPTION: Wavefront Alignment Algorithms benchmarking tool
 */

#include <omp.h>

#include "align_benchmark_params.h"

#include "utils/commons.h"
#include "utils/sequence_buffer.h"
#include "system/profiler_timer.h"

#include "alignment/score_matrix.h"
#include "edit/edit_dp.h"
#include "gap_linear/nw.h"
#include "gap_affine/swg.h"
#include "gap_affine2p/affine2p_matrix.h"
#include "gap_affine2p/affine2p_dp.h"
#include "wavefront/wavefront_align.h"

#include "benchmark/benchmark_indel.h"
#include "benchmark/benchmark_edit.h"
#include "benchmark/benchmark_gap_linear.h"
#include "benchmark/benchmark_gap_affine.h"
#include "benchmark/benchmark_gap_affine2p.h"
#include "benchmark/benchmark_check.h"

#include "../foralign/hirschberg.h"

/*
 * WFA lambda (custom match function)
 */
typedef struct {
  char* pattern;
  int pattern_length;
  char* text;
  int text_length;
} match_function_params_t;
match_function_params_t lambda_params;
// Simplest Extend-matching function (for testing purposes)
int lambda_function(int v,int h,void* arguments) {
  // Extract parameters
  match_function_params_t* const match_arguments = (match_function_params_t*)arguments;
  // Check match
  if (v >= match_arguments->pattern_length || h >= match_arguments->text_length) return 0;
  return (match_arguments->pattern[v] == match_arguments->text[h]);
}
/*
 * Algorithms
 */
bool align_benchmark_is_wavefront(
    const alignment_algorithm_type algorithm) {
  return algorithm == alignment_indel_wavefront ||
         algorithm == alignment_edit_wavefront ||
         algorithm == alignment_gap_linear_wavefront ||
         algorithm == alignment_gap_affine_wavefront ||
         algorithm == alignment_gap_affine2p_wavefront;
}
/*
 * Configuration
 */
wavefront_aligner_t* align_input_configure_wavefront(
    align_input_t* const align_input) {
  // Set attributes
  wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
  attributes.memory_mode = parameters.wfa_memory_mode;
  if (parameters.wfa_score_only) {
    attributes.alignment_scope = compute_score;
  }
  // WF-Heuristic
  switch (parameters.wfa_heuristic) {
    case wf_heuristic_none:
      attributes.heuristic.strategy = wf_heuristic_none;
      break;
    case wf_heuristic_banded_static:
      attributes.heuristic.strategy = wf_heuristic_banded_static;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_banded_adaptive:
      attributes.heuristic.strategy = wf_heuristic_banded_adaptive;
      attributes.heuristic.min_k = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_k = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_wfadaptive:
      attributes.heuristic.strategy = wf_heuristic_wfadaptive;
      attributes.heuristic.min_wavefront_length = parameters.wfa_heuristic_p1;
      attributes.heuristic.max_distance_threshold = parameters.wfa_heuristic_p2;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p3;
      break;
    case wf_heuristic_xdrop:
      attributes.heuristic.strategy = wf_heuristic_xdrop;
      attributes.heuristic.xdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    case wf_heuristic_zdrop:
      attributes.heuristic.strategy = wf_heuristic_zdrop;
      attributes.heuristic.zdrop = parameters.wfa_heuristic_p1;
      attributes.heuristic.steps_between_cutoffs = parameters.wfa_heuristic_p2;
      break;
    default:
      break;
  }
  // Select flavor
  switch (parameters.algorithm) {
    case alignment_indel_wavefront:
      attributes.distance_metric = indel;
      break;
    case alignment_edit_wavefront:
      attributes.distance_metric = edit;
      break;
    case alignment_gap_linear_wavefront:
      attributes.distance_metric = gap_linear;
      attributes.linear_penalties = parameters.linear_penalties;
      break;
    case alignment_gap_affine_wavefront:
      attributes.distance_metric = gap_affine;
      attributes.affine_penalties = parameters.affine_penalties;
      break;
    case alignment_gap_affine2p_wavefront:
      attributes.distance_metric = gap_affine_2p;
      attributes.affine2p_penalties = parameters.affine2p_penalties;
      break;
    default:
      return NULL; // No WF selected
      break;
  }
  // Select alignment form
  if (parameters.align_span_extension) {
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.extension = true;
  } else if (parameters.align_span_endsfree) {
    attributes.alignment_form.span = alignment_endsfree;
    attributes.alignment_form.extension = false;
  } else { // Global
    attributes.alignment_form.span = alignment_end2end;
    attributes.alignment_form.extension = false;
  }
  // Misc
  attributes.plot.enabled = (parameters.plot != 0);
  attributes.plot.align_level = (parameters.plot < 0) ? -1 : parameters.plot - 1;
  attributes.system.verbose = parameters.verbose;
  attributes.system.max_memory_abort = parameters.wfa_max_memory;
  attributes.system.max_alignment_steps = parameters.wfa_max_steps;
  attributes.system.max_num_threads = parameters.wfa_max_threads;
  // Allocate
  return wavefront_aligner_new(&attributes);
}
void align_input_configure_global(
    align_input_t* const align_input) {
  // Clear
  benchmark_align_input_clear(align_input);
  // Penalties
  align_input->linear_penalties = parameters.linear_penalties;
  align_input->affine_penalties = parameters.affine_penalties;
  align_input->affine2p_penalties = parameters.affine2p_penalties;
  // ref input
  align_input->input_ref_file = parameters.input_ref_file;
  // Output
  align_input->output_file = parameters.output_file;
  align_input->output_full = parameters.output_full;
  // MM
  align_input->mm_allocator = mm_allocator_new(BUFFER_SIZE_1M);
  // WFA
  if (align_benchmark_is_wavefront(parameters.algorithm)) {
    if (parameters.wfa_lambda) {
      align_input->wfa_match_funct = lambda_function;
      align_input->wfa_match_funct_arguments = &lambda_params;
    }
    align_input->wf_aligner = align_input_configure_wavefront(align_input);
  } else {
    align_input->wf_aligner = NULL;
  }
  // PROFILE/STATS
  timer_reset(&align_input->timer);
  // DEBUG
  align_input->debug_flags = 0;
  align_input->debug_flags |= parameters.check_metric;
  if (parameters.check_display) align_input->debug_flags |= ALIGN_DEBUG_DISPLAY_INFO;
  if (parameters.check_correct) align_input->debug_flags |= ALIGN_DEBUG_CHECK_CORRECT;
  if (parameters.check_score) align_input->debug_flags |= ALIGN_DEBUG_CHECK_SCORE;
  if (parameters.check_alignments) align_input->debug_flags |= ALIGN_DEBUG_CHECK_ALIGNMENT;
  align_input->check_linear_penalties = &parameters.linear_penalties;
  align_input->check_affine_penalties = &parameters.affine_penalties;
  align_input->check_affine2p_penalties = &parameters.affine2p_penalties;
  align_input->check_bandwidth = parameters.check_bandwidth;
  align_input->verbose = parameters.verbose;
}
void align_input_configure_local(
    align_input_t* const align_input) {
  // Ends-free configuration
  if (parameters.align_span_endsfree) {
    const int plen = align_input->pattern_length;
    const int tlen = align_input->text_length;
    align_input->pattern_begin_free = nominal_prop_u32(plen,parameters.pattern_begin_free);
    align_input->pattern_end_free = nominal_prop_u32(plen,parameters.pattern_end_free);
    align_input->text_begin_free = nominal_prop_u32(tlen,parameters.text_begin_free);
    align_input->text_end_free = nominal_prop_u32(tlen,parameters.text_end_free);
    if (align_benchmark_is_wavefront(parameters.algorithm)) {
      wavefront_aligner_set_alignment_free_ends(align_input->wf_aligner,
          align_input->pattern_begin_free,align_input->pattern_end_free,
          align_input->text_begin_free,align_input->text_end_free);
    }
  }
  // Custom extend-match function
  if (parameters.wfa_lambda) {
    lambda_params.pattern = align_input->pattern;
    lambda_params.pattern_length = align_input->pattern_length;
    lambda_params.text = align_input->text;
    lambda_params.text_length = align_input->text_length;
  }
  align_input->debug_flags = true;
}
void align_benchmark_free(
    align_input_t* const align_input) {
  if (align_input->wf_aligner) wavefront_aligner_delete(align_input->wf_aligner);
  mm_allocator_delete(align_input->mm_allocator);
}
/*
 * I/O
 */
bool align_benchmark_read_ref_input(
    FILE* input_file,
    char** line,
    size_t* line_allocated,
    align_input_t* const align_input) {
  // Parameters
  int line_length=0;
  // Read > line
  line_length = getline(line,line_allocated,input_file);
  if (line_length==-1) return false;
  // Read sequence line
  line_length = getline(line,line_allocated,input_file);
  if (line_length==-1) return false;
  // Configure input
  align_input->sequence_id = 1;
  align_input->pattern = *line;
  align_input->pattern_length = line_length - 1;
  if (align_input->pattern[align_input->pattern_length] == '\n') {
    align_input->pattern[align_input->pattern_length] = '\0';
  }
  return true;
}

bool align_benchmark_read_seq_input(
    FILE* input_file,
    char** line,
    size_t* line_allocated,
    align_input_t* const align_input) {
  // Parameters
  int line_length=0;
  // Read > line
  line_length = getline(line,line_allocated,input_file);
  if (line_length==-1) return false;
  // Read sequence line
  line_length = getline(line,line_allocated,input_file);
  if (line_length==-1) return false;
  // Configure input
  align_input->text = *line;
  align_input->text_length = line_length - 1;
  if (align_input->text[align_input->text_length] == '\n') {
    align_input->text[align_input->text_length] = '\0';
  }
  return true;
}


/*
 * Display
 */

void print_Error(
    FILE* const stream,
    align_input_t* const align_input,
    const bool print_wf_stats) {
  fprintf(stream,"Aligned bases rate: \n %.2f %%\n", 100.0 * align_input->align_matches.total / align_input->align_bases.total);
  // CIGAR stats
  fprintf(stream," Details:                  \n");
  fprintf(stream,"   => CIGAR.Matches        ");
  counter_print(stream,&align_input->align_matches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Mismatches     ");
  counter_print(stream,&align_input->align_mismatches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Insertions     ");
  counter_print(stream,&align_input->align_ins,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Deletions      ");
  counter_print(stream,&align_input->align_del,&align_input->align_bases,"bases     ",true);
}


/*
 * Benchmark
 */
void align_benchmark_run_algorithm(
    align_input_t* const align_input) {
  // Sequence-dependent configuration
  align_input_configure_local(align_input);
  // Select algorithm
  switch (parameters.algorithm) {
    // Indel
    case alignment_indel_wavefront:
      benchmark_indel_wavefront(align_input);
      break;
    // Edit
    case alignment_edit_bpm:
      benchmark_edit_bpm(align_input);
      break;
    case alignment_edit_dp:
      benchmark_edit_dp(align_input);
      break;
    case alignment_edit_dp_banded:
      benchmark_edit_dp_banded(align_input,parameters.bandwidth);
      break;
    case alignment_edit_wavefront:
      benchmark_edit_wavefront(align_input);
      break;
    // Gap-linear
    case alignment_gap_linear_nw:
      benchmark_gap_linear_nw(align_input,&parameters.linear_penalties);
      break;
    case alignment_gap_linear_wavefront:
      benchmark_gap_linear_wavefront(align_input,&parameters.linear_penalties);
      break;
    // Gap-affine
    case alignment_gap_affine_swg:
      benchmark_gap_affine_swg(align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_endsfree:
      benchmark_gap_affine_swg_endsfree(
          align_input,&parameters.affine_penalties);
      break;
    case alignment_gap_affine_swg_banded:
      benchmark_gap_affine_swg_banded(align_input,
          &parameters.affine_penalties,parameters.bandwidth);
      break;
    case alignment_gap_affine_wavefront:
      benchmark_gap_affine_wavefront(align_input,&parameters.affine_penalties);
      break;
    // Gap-affine 2p
    case alignment_gap_affine2p_dp:
      benchmark_gap_affine2p_dp(align_input,&parameters.affine2p_penalties);
      break;
    case alignment_gap_affine2p_wavefront:
      benchmark_gap_affine2p_wavefront(align_input,&parameters.affine2p_penalties);
      break;
    case alignment_hirschberg_multi:
    {
      cigar_t* const cigar = cigar_new(align_input->pattern_length+align_input->text_length);
      // Align
      timer_start(&align_input->timer);
      hirschberg_cigar(align_input->pattern, align_input->text, align_input->pattern_length, align_input->text_length,
                       parameters.affine_penalties.gap_opening, parameters.affine_penalties.gap_extension,
                       parameters.affine_penalties.match, parameters.affine_penalties.mismatch, 0, DNA,
                       parameters.divide_threads, cigar->operations, &(cigar->end_offset), 1, 1, parameters.dp_threads);
      cigar->begin_offset = 0;
      timer_stop(&align_input->timer);
      
      counter_add(&(align_input->align_bases),align_input->pattern_length);
      int i;
      for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
        switch (cigar->operations[i]) {
          case 'M': counter_add(&(align_input->align_matches),1); break;
          case 'X': counter_add(&(align_input->align_mismatches),1); break;
          case 'I': counter_add(&(align_input->align_ins),1); break;
          case 'D': default: counter_add(&(align_input->align_del),1); break;
        }
      }

      cigar_free(cigar);
    }
      break;
    case alignment_russians:
    {
      cigar_t* const cigar = cigar_new(align_input->pattern_length+align_input->text_length);
      // Align
      timer_start(&align_input->timer);
      hirschberg_cigar(align_input->pattern, align_input->text, align_input->pattern_length, align_input->text_length,
                       parameters.affine_penalties.gap_opening, parameters.affine_penalties.gap_extension,
                       parameters.affine_penalties.match, parameters.affine_penalties.mismatch, parameters.russian_block_size,
                       DNA, parameters.divide_threads, cigar->operations, &(cigar->end_offset), 0, 1, parameters.dp_threads);
      cigar->begin_offset = 0;
      timer_stop(&align_input->timer);
      
      counter_add(&(align_input->align_bases),align_input->pattern_length);
      int i;
      for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
        switch (cigar->operations[i]) {
          case 'M': counter_add(&(align_input->align_matches),1); break;
          case 'X': counter_add(&(align_input->align_mismatches),1); break;
          case 'I': counter_add(&(align_input->align_ins),1); break;
          case 'D': default: counter_add(&(align_input->align_del),1); break;
        }
      }

      cigar_free(cigar);
    }
      break;
    default:
      fprintf(stderr,"Algorithm not implemented\n");
      exit(1);
      break;
  }
}
void align_benchmark_sequential() {
  // PROFILE
  timer_reset(&parameters.timer_global);
  // I/O files
  parameters.input_file = fopen(parameters.input_filename, "r");
  parameters.input_ref_file = fopen(parameters.input_ref_filename, "r");
  if (parameters.input_file == NULL || parameters.input_ref_file == NULL) {
    fprintf(stderr,"Input file '%s' or '%s' couldn't be opened\n",parameters.input_filename, parameters.input_ref_filename);
    exit(1);
  }
  // Global configuration
  align_input_t align_input;
  align_input_configure_global(&align_input);
  // Read-align loop
  // Read input sequence-pair
  bool input_read = align_benchmark_read_ref_input(parameters.input_ref_file,&parameters.line1,&parameters.line1_allocated,&align_input);
  input_read &=  align_benchmark_read_seq_input(parameters.input_file,&parameters.line2,&parameters.line2_allocated,&align_input);
  if (!input_read) { fprintf(stderr, "Warning: no sequence found in file. Will generate inaccurate values.\n"); }
  // Execute the selected algorithm
  timer_start(&parameters.timer_global); // PROFILE
  align_benchmark_run_algorithm(&align_input);
  timer_stop(&parameters.timer_global); // PROFILE
  print_Error(stderr,&align_input,true);
  // Free
  align_benchmark_free(&align_input);
  fclose(parameters.input_file);
  if (parameters.output_file) fclose(parameters.output_file);
  free(parameters.line1);
  free(parameters.line2);
}
/*
 * Main
 */
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  align_benchmark_sequential();
  return 0;
}
