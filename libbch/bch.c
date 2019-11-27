#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "gf2x.h"

#include "additive_fft.h"
#include "multiplicative_fft.h"
#include "additive_fft.h"
#include "finite_field.h"
#include "gf2_polynomials.h"
#include "gf2_polynomials_misc.h"

#include "bch.h"

static const uint32_t sUL = sizeof(unsigned long);

#define LOGGING
#define USE_ADDITIVE_DFT
#define USE_ADDITIVE_DFT_SYN

void multiply_polynomial_by_monomial(
    TfiniteField* p_f,
    uint32_t* p_poly, uint32_t p_polyMaxDegree,
    uint32_t p_rootLog);

TbchCode bch_gen_code(uint32_t p_codeword_size, uint32_t p_minimum_distance)
{
  TbchCode code;
  code.m_distance = 0;
  code.m_packed_gen_poly = NULL;
  code.m_length = 0;
  code.m_gen_poly_degree = 0;
  if((p_minimum_distance >= p_codeword_size) || p_codeword_size > (1u << g_maxDegree)) return code;

  uint32_t degree = 2;
  while ((1u << degree) < p_codeword_size) degree++;

  TfiniteField f = create_binary_finite_field(degree);
  code = bch_build_polynomial(&f, p_minimum_distance);
  if (code.m_gen_poly_degree >= p_codeword_size)
  {
    // error: the code has no data bits, because there are too many redundancy bits
    bch_free_code(&code);
  }
  else
  {
    code.m_length = p_codeword_size;
  }
  return code;
}

void bch_free_code(TbchCode *p_code)
{
  p_code->m_distance = 0;
  free_finite_field(&(p_code->m_field));
  free(p_code->m_packed_gen_poly);
  p_code->m_packed_gen_poly = NULL;
}

TbchCode bch_build_polynomial(TfiniteField* p_f, uint32_t p_minimumDistance)
{
  TbchCode code;
  code.m_field = *p_f;
  code.m_distance = p_minimumDistance;
  code.m_length = 0; // to be set later, to a value greater than the number of redundancy bits
  // (the degree of the generating polynomial, m_degree)

  // the array below stores the minimal polynomial of each x**i with i in the code interval,
  // in compact form, in its first half. Second half is to store the product of these polynomials
  unsigned long* polynomials = (unsigned long*) calloc(2 * p_minimumDistance, sUL);

#ifdef LOGGING
  printf("finding redundancies in zeros of generating polynomial...\n");
#endif
  uint32_t numPolys;
  // 'remaining' keeps track of the values belonging to the cycles already discovered
  uint32_t* remaining = (uint32_t*) malloc(p_f->n * 4);
  bch_explore_cycles(
        p_f, p_minimumDistance, 0, remaining, polynomials, &numPolys,
        &code.m_gen_poly_degree);
  free(remaining);
#ifdef LOGGING
  printf("Code has %d redundancy bits\n", code.m_gen_poly_degree);
  printf("...done\n");
#endif
  uint32_t dir = 0;
  uint32_t stride = 1;
#ifdef LOGGING
  printf("Multiplying the minimal polynomials to get the generating polynomial...\n");
#endif
  while(numPolys > 1)
  {
    unsigned long* source = polynomials + dir       * p_minimumDistance;
    unsigned long* dest   = polynomials + (1 - dir) * p_minimumDistance;
    for(uint32_t i = 0; i < numPolys / 2; i++)
    {
      gf2x_mul(
            dest + 2 * i * stride,
            source + 2 * i *stride,
            stride,
            source + (2 * i + 1) * stride,
            stride
            );
    }
    if(numPolys & 1) memcpy(dest + (numPolys - 1)*stride, source + (numPolys - 1)*stride, stride *sUL);
    numPolys = (numPolys + 1) / 2;
    stride = stride * 2;
    dir = 1 - dir;
  }
  const uint32_t packedULsize = degree_to_packed_size_u32(code.m_gen_poly_degree);
  code.m_packed_gen_poly = (unsigned long*) calloc(packedULsize, sUL);
  memcpy(code.m_packed_gen_poly, polynomials + dir * p_minimumDistance, packedULsize * sUL);
  free(polynomials);
  return code;
}

void bch_find_best_offset(TfiniteField* p_f, uint32_t p_minimum_distance)
{
  uint32_t *remaining;
  uint32_t offset, polyDegree = 0, minPolyDegree = 0, minArgOffset;

  minArgOffset = 0;
  remaining = (uint32_t*) malloc(p_f->n * sizeof (uint32_t));
  for (offset = 0; offset < p_f->n - 1 - p_minimum_distance; offset++)
  {
    bch_explore_cycles(p_f, p_minimum_distance, offset, remaining, 0, 0, &polyDegree);
    printf("Offset %d, degree %d\n", offset, polyDegree);
    if (offset == 0 || minPolyDegree > polyDegree)
    {
      minPolyDegree = polyDegree;
      minArgOffset = offset;
    }
  }
  printf("Degree with offset %d: %d\n", minArgOffset, minPolyDegree);
  free(remaining);
}


void bch_explore_cycles(TfiniteField* p_f, uint32_t p_numValues, uint32_t p_offset,
    uint32_t* p_remainingPtr, unsigned long* p_binaryPolynomials, uint32_t* po_numPolys,
    uint32_t* po_totalDegree)
{
  uint32_t elementDegree, cycleTagged, genPolyDegree, elementLog, initialElementLog;
  uint32_t numFoundValues, i, sz;

  sz = p_f->n - 1;
  numFoundValues = 0;
  //mark 0
  p_remainingPtr[0] = 1;

  uint32_t cyclePoly[g_maxDegree + 1];
  uint32_t polynomialIndex = 0;

  for (i = 1; i < sz + 1; i++)
  {
    // invariant on p_remainingPtr:
    // if p_remainingPtr[i] = i then i has not yet been visited as part of a cycle exploration
    // otherwise p_remainingPtr[i] > i and by iterating i = p_remainingPtr[i],
    // and stopping at a fixed point, one finds the first i s.t. i has not yet been visited
    // or one exits (when i = p_sz)
    p_remainingPtr[i] = i;
  }
  initialElementLog = 1;
  genPolyDegree = 0;
  do
  {
    // walk through a new cycle
    elementDegree = 0;
    cycleTagged = 0;
    elementLog = initialElementLog;
    cyclePoly[0] = 1;
    do
    {
      if ((elementLog >= p_offset) && (elementLog < p_numValues + p_offset))
      {
        // a root was found in the current cycle
        numFoundValues++;
        cycleTagged = 1;
      }
      multiply_polynomial_by_monomial(p_f, cyclePoly, elementDegree++, elementLog);
      // mark element to avoid considering it later
      p_remainingPtr[elementLog] = p_remainingPtr[elementLog + 1];
      elementLog = (2 * elementLog) % sz;
    }
    while (elementLog != initialElementLog);

    if (cycleTagged)
    {
      //if(elementDegree < p_finiteFieldPtr->m) printf("%d ", elementDegree);
      if (p_binaryPolynomials)
      {
        p_binaryPolynomials[polynomialIndex] = 0;
        //degree = 1 + num. of roots
        for (i = 0; i < elementDegree + 1; i++)
        {
          set_bit(&p_binaryPolynomials[polynomialIndex], i, cyclePoly[i]);
        }
        polynomialIndex++;
      }
      genPolyDegree += elementDegree;
    }
    // Find next cycle initial point (first point not encountered in previous cycles)
    while (initialElementLog != sz && initialElementLog < p_remainingPtr[initialElementLog])
    {
      initialElementLog = p_remainingPtr[initialElementLog];
      // b = x**initialElementLog is a new element (not seen in previous cycles)
      // hence its conjugates are also new
      // (otherwise a conjugate a of b would have already been obtained, but then b would be a
      // conjugate of a and would have been already been found)
      // this is why we don't test whether the elements in the cycle of b are new
    }
  }
  while (initialElementLog != sz && numFoundValues < p_numValues);
  if(po_totalDegree) *po_totalDegree = genPolyDegree;
  if(po_numPolys) *po_numPolys = polynomialIndex;
}

void bch_reduce_syndrome(TbchCode *p_code, unsigned long* pio_encoded_data)
{
  euclidean_division_packed_inplace_gf2(pio_encoded_data, p_code->m_length - 1,
      p_code->m_packed_gen_poly, p_code->m_gen_poly_degree, 0);
}

int32_t bch_eval_syndrome(TbchCode* p_code, TbchCodeBuffers* p_buffer, unsigned long* p_packed_syndrome)
{
  uint32_t i, j;
  int32_t syn_error = 0;
  TfiniteField* f = &p_code->m_field;
  const uint32_t* exp = f->exp;
  const uint32_t order = f->n - 1;
  uint32_t* syndrome = p_buffer->m_syndrome;
  uint32_t* buffer = p_buffer->m_buffer_1;

  //2. constant : empirically chosen to equalize Chien search and FFT running times at threshold
#ifndef USE_ADDITIVE_DFT_SYN
  float FFT_threshold = (float) 2. * getFactorSum(f);
#else
  float l = logf((float) f->m);
  float FFT_threshold = 1.f * f->m * l * l;
#endif
  float cost =
      ((float) p_code->m_distance) *
      ((float) p_code->m_gen_poly_degree) / ((float) p_code->m_length);
#ifdef LOGGING
  printf("  Syndrome computation: must evaluate polynomial of degree %d at %d points; "
         "cost: %.0f; FFT threshold: %.0f\n", p_code->m_gen_poly_degree,
         p_code->m_distance, (double) cost, (double) FFT_threshold);
#endif

  int use_fft = cost > FFT_threshold;
  int pre_reduce_syndrome = !use_fft;

  uint32_t max_idx = 0;
  unsigned long* source_buffer = NULL;
  if(pre_reduce_syndrome)
  {
    // reduce codeword by syndrome polynomial first
    memcpy(buffer, p_packed_syndrome, degree_to_packed_size(p_code->m_length - 1) * sUL);
    bch_reduce_syndrome(p_code, (unsigned long*) p_buffer->m_buffer_1);
    max_idx = p_code->m_gen_poly_degree;
    source_buffer = (unsigned long*) buffer;
  }
  else
  {
    max_idx = p_code->m_length;
    source_buffer = p_packed_syndrome;
  }

  if (use_fft)
  {
#ifdef LOGGING
    printf("  Using FFT to evaluate syndrome polynomial\n");
#endif
    printf("max_idx before reduction: %d\n", max_idx);
    for (i = 0; i < max_idx; i++) syndrome[i] = get_bit(source_buffer, i);
    for (     ; i < f->n;    i++) syndrome[i] = 0;

#ifndef USE_ADDITIVE_DFT_SYN
    if(max_idx == f->n)
    {
      //multiplicative FFT only takes into account terms of degree 0 ... order - 1: fold term of degree 'order'
      syndrome[0] ^= syndrome[order];
      syndrome[order] = 0;
      max_idx = order;
    }
    printf("max_idx: %d\n", max_idx);
    evaluate_polynomial_multiplicative_FFT(f, syndrome, syndrome, max_idx, buffer);
#else
    evaluate_polynomial_additive_FFT(f, syndrome, buffer, max_idx);
    for(i = 0; i < f->n; i++) syndrome[i] = buffer[i];
#endif
  }
  else
  {
#ifdef LOGGING
    printf("  Using Chien search to evaluate syndrome polynomial\n");
#endif
    for (i = 0; i < p_code->m_distance; i++) syndrome[i] = 0;
    for (j = 0; j < max_idx; j++)
    {
      if (get_bit(source_buffer, j))
      {
        uint32_t ij = 0;
        for (i = 0; i < p_code->m_distance; i++)
        {
          syndrome[i] ^= exp[ij];
          ij += j;
          if (ij >= order) ij -= order;
        }
      }
    }
  }

  for (i = 0; i < p_code->m_distance; i++)
  {
    if (syndrome[i] != 0)
    {
      syn_error = 1;
      break;
    }
  }
  return syn_error;
}

uint32_t bch_compute_elp(TbchCode* p_code, TbchCodeBuffers* p_buffer)
{
  TfiniteField* f = &p_code->m_field;
  const uint32_t order = f->n - 1;
  uint32_t* syndrome = p_buffer->m_syndrome;
  uint32_t* A = p_buffer->m_A;
  uint32_t* B = p_buffer->m_B;
  uint32_t* C = p_buffer->m_C_elp;

  uint32_t Adegree = 0, Cdegree = 0, Bdegree = 0;
  B[0] = 1;
  C[0] = 1;
  uint32_t inv_b_log = 0; // initially b = 1
  uint32_t L = 0;
  uint32_t m = 1;
  for (uint32_t i = 0; i < p_code->m_distance - 1; i++)
  {
    uint32_t j;
    // calculate discrepancy d_{i+1}
    uint32_t d = syndrome[i + 1];
    for (j = 1; j < Cdegree + 1; j++) d ^= multiply(f, C[j], syndrome[i + 1 - j]);
    if (d == 0) m = m + 1; // annihilation continues
    else
    {
      // A(x) = C(x) - (d / b) x**m . B(x)
      // compute log of leading term of B(x)
      uint32_t d_log = f->log[d];
      uint32_t mult_log = d_log + inv_b_log;
      if(mult_log >= order) mult_log -= order;
      // compute A and its degree
      uint32_t min_m_cterms = m < Cdegree + 1 ? m : Cdegree + 1; // min_m_cterms <= m
      for(j = 0; j < min_m_cterms   ; j++) A[j] = C[j];
      for(     ; j < m              ; j++) A[j] = 0;
      for(     ; j < Bdegree + m + 1; j++) A[j] = multiply_by_log(f, B[j - m], mult_log);
      for(     ; j < Cdegree + 1    ; j++) A[j] = 0;
      if(Cdegree != Bdegree + m)
      {
        Adegree = Cdegree > Bdegree + m ? Cdegree : Bdegree + m;
        for (j = m; j < Cdegree + 1; j++) A[j] ^= C[j];
      }
      else
      {
        for(j = m; j < Cdegree + 1; j++)
        {
          A[j] ^= C[j];
          if(A[j]) Adegree = j;
        }
      }

      uint32_t* swap = NULL;
      if(2 * L <= i)
      {
        //swap B and C
        // (we are only interested in B = C, but need to track other pointers for buffer reuse)
        L = i + 1 - L;
        Bdegree = Cdegree;
        swap = B; B = C; C = swap;
        inv_b_log = order - d_log;
        m = 1;
      }
      else
      {
        m = m + 1;
      }
      //swap A and C
      // (we are only interested in C = A, but need to track other pointers for buffer reuse)
      swap = C; C = A; A = swap;
      Cdegree = Adegree;
    }
  }
  // swap buffers in order to ensure that the ELP is in p_buffer->m_elp
  // this is possible since all buffers have the same size
  p_buffer->m_C_elp = C;
  p_buffer->m_A = A;
  p_buffer->m_B = B;
  return L;
}

uint32_t bch_locate_errors(TbchCode* p_code, TbchCodeBuffers* p_buffer, uint32_t p_elp_degree)
{
  TfiniteField* f = &p_code->m_field;
  uint32_t order = f->n - 1;
  const uint32_t* log_t = f->log;
  const uint32_t* exp_t = f->exp;
  uint32_t* syndrome = p_buffer->m_syndrome;
  uint32_t* result_buffer = p_buffer->m_buffer_1;
  uint32_t* aux_buffer = p_buffer->m_buffer_2;
  uint32_t* elp = p_buffer->m_C_elp;

  //check if the degree of the elp is big enough to use the fft
#ifdef LOGGING
  printf("  Degree of error locator polynomial: %d\n", p_elp_degree);
#endif
  // constants below empirically chosen to equalize Chien search and FFT running times at threshold
  // the constants are different than in the syndrome computation case because here the polynomial
  // is in GF(2^m), whereas for the syndrome it is on GF(2). This does not change the FFT
  // computation time (always done in GF(2^m)), but does change the Chien search computation time.

#ifndef USE_ADDITIVE_DFT
  float FFT_threshold = 0.8f * (float) (getFactorSum(f));
#else
  float l = logf((float) f->m);
  float FFT_threshold = 0.5f * f->m * l * l;
#endif
#ifdef LOGGING
  printf("  Error location: has degree %d, FFT threshold: %.0f\n", p_elp_degree, (double) FFT_threshold);
#endif

  int use_fft = p_elp_degree > FFT_threshold;

  if (use_fft)
  {
#ifdef LOGGING
    printf("  Using FFT to find roots of ELP\n");
#endif
    {
      uint32_t i;
      for (i = 0; i < p_elp_degree + 1; i++) syndrome[i] = elp[i];
      for (     ; i < order           ; i++) syndrome[i] = 0;
    }
#ifndef USE_ADDITIVE_DFT
    evaluate_polynomial_multiplicative_FFT(f, syndrome, result_buffer, p_elp_degree + 1, aux_buffer);
#else
    evaluate_polynomial_additive_FFT(f, syndrome, result_buffer, p_elp_degree + 1);
#endif
  }
  else
  {
#ifdef LOGGING
    printf("  Using Chien search to find roots of ELP\n");
#endif
    // put elp into logarithm form
    for (uint32_t j = 0; j < p_elp_degree + 1; j++) aux_buffer[j] = log_t[elp[j]];
    for (uint32_t i = 0; i < order; i++)
    {
      // evaluate ELP on x**i
      uint32_t v = elp[0];
      for (uint32_t j = 1; j < p_elp_degree + 1; j++)
      {
        uint32_t log_coeff = aux_buffer[j];
        if (log_coeff != f->logzero)
        {
          v ^= exp_t[log_coeff];
          log_coeff += j;
          if(log_coeff >= order) log_coeff -= order;
          aux_buffer[j] = log_coeff;
        }
      }
      //ELP(x**i)
      result_buffer[i] = v;
    }
  }
  uint32_t count = 0;
  for (uint32_t i = 0; i < order; i++)
  {
    if (result_buffer[i] == 0)
    {
      uint32_t errorIdx = order - i;
      //if(errorIdx < p_codePtr->m_length)
      p_buffer->m_error[count++] = errorIdx;
    }
  }
  return count;
}

int bch_encode(TbchCode *p_code, unsigned long* po_codeword, uint32_t* p_data, uint32_t p_data_size)
{
  const uint32_t d = p_code->m_gen_poly_degree;
  if(p_data_size > p_code->m_length - d) return 1;
  memset(po_codeword, 0, degree_to_packed_size(p_code->m_length - 1) * sUL);
  // copy X**(code.m_gen_poly_degree) * 'data',  where 'data' is seen as
  // a polynomial, into encodedData...
  for (uint32_t i = 0; i < p_data_size; i++) set_bit(po_codeword, i + d, p_data[i]);
  // ... reduce modulo the generating polynomial ...
  bch_reduce_syndrome(p_code, po_codeword);
  // then form the full codeword by copying data again.
  // the resulting polynomial cancels on all roots of the generating polynomial.
  for (uint32_t i = 0; i < p_data_size; i++) set_bit(po_codeword, i + d, p_data[i]);
  return 0;
}

uint32_t bch_decode(TbchCode* p_code, unsigned long* pio_packed_data,
    TbchCodeBuffers* p_buffer, uint32_t* po_corrected)
{
  uint32_t num_errors_found = 0;
  uint32_t ok = 0;
  TbchCodeBuffers buffer_struct;

  if (p_buffer) buffer_struct = *p_buffer;
  else          buffer_struct = bch_alloc_buffers(p_code);

#ifdef LOGGING
  printf("Evaluating syndrome on roots of the generator polynomial...\n");
#endif
  int32_t syn_error = bch_eval_syndrome(p_code, &buffer_struct, pio_packed_data);
#ifdef LOGGING
  printf("... done\n");
#endif
  if (syn_error)
  {
#ifdef LOGGING
    printf("Computing error location polynomial...\n");
#endif
    uint32_t deg = bch_compute_elp(p_code, &buffer_struct);
    if (deg <= (p_code->m_distance - 1) / 2)
    {
      // Find the roots of the ELP that are in the index range
      // 0 ... p_codePtr->m_length - 1
      num_errors_found = bch_locate_errors(p_code, &buffer_struct, deg);
      if (num_errors_found == deg)
      {
        ok = 1;
        for (uint32_t i = 0; i < num_errors_found; i++)
        {
          add_bit(pio_packed_data, buffer_struct.m_error[i]);
        }
      }
#ifdef LOGGING
      else
      {
        printf("BCH: elp has degree %d, but found only %d roots\n", deg, num_errors_found);
      }
#endif
    }
#ifdef LOGGING
    else
    {
      printf("BCH: Too many errors, can't correct\n");
    }
#endif
  }
  else ok = 1;
  if (!p_buffer) bch_free_buffers(&buffer_struct);
  if (po_corrected) *po_corrected = num_errors_found;
  return ok;
}

TbchCodeBuffers bch_alloc_buffers(TbchCode* p_code)
{
  TbchCodeBuffers buffer;

  uint32_t byteSize = (p_code->m_distance + 1) * 4;
  buffer.m_syndrome = (uint32_t*) malloc(p_code->m_field.n * 4);
  buffer.m_buffer_1 = (uint32_t*) malloc(p_code->m_field.n * 4);
  buffer.m_buffer_2 = (uint32_t*) malloc(p_code->m_field.n * 4);
  buffer.m_A = (uint32_t*) malloc(byteSize);
  buffer.m_B = (uint32_t*) malloc(byteSize);
  buffer.m_C_elp = (uint32_t*) malloc(byteSize);
  buffer.m_error = (uint32_t*) malloc(byteSize);
  return buffer;
}

void bch_free_buffers(TbchCodeBuffers* p_buffer)
{
  free(p_buffer->m_C_elp);
  p_buffer->m_C_elp = 0;
  free(p_buffer->m_A);
  p_buffer->m_A = 0;
  free(p_buffer->m_B);
  p_buffer->m_B = 0;
  free(p_buffer->m_syndrome);
  p_buffer->m_syndrome = 0;
  free(p_buffer->m_error);
  p_buffer->m_error = 0;
  free(p_buffer->m_buffer_2);
  p_buffer->m_buffer_2 = 0;
  free(p_buffer->m_buffer_1);
  p_buffer->m_buffer_1 = 0;
}

