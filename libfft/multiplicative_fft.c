#include "multiplicative_fft.h"
#include "finite_field.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//factors of 2^8-1 to 2^26-1
static const uint32_t minLogSize = 8;
static const uint32_t maxLogSize = 26;
static uint32_t fft_factors[] = {
  3, 5, 17,                    // 2^8 - 1
  7, 73,                       // 2^9 - 1
  3, 11, 31,                   // 2^10 - 1
  23, 89,                      // 2^11 - 1
  3, 3, 5, 7, 13,              // 2^12 - 1
  8191,                        // 2^13 - 1
  3, 43, 127,                  // 2^14 - 1
  7, 31, 151,                  // 2^15 - 1
  3, 5, 17, 257,               // 2^16 - 1
  131071,                      // 2^17 - 1
  3, 3, 3, 7, 19, 73,          // 2^18 - 1
  524287,                      // 2^19 - 1
  3, 5, 5, 11, 31, 41,         // 2^20 - 1
  7, 7, 127, 337,              // 2^21 - 1
  3, 23, 89, 683,              // 2^22 - 1
  47, 178481,                  // 2^23 - 1
  3, 3, 5, 7, 13, 17, 241,     // 2^24 - 1
  31, 601, 1801,               // 2^25 - 1
  3, 2731, 8191};              // 2^26 - 1

static uint32_t fft_offset[] = {
  0, 3, 5, 8, 10, 15, 16, 19, 22, 26, 27, 33, 34, 40, 44, 48, 50, 57, 60, 63};

uint32_t factorOffset(uint32_t p_log)
{
  return fft_offset[p_log - minLogSize];
}

uint32_t getFactorSum(TfiniteField* p_f)
{
  uint32_t l_numFactors, l_sum, i;

  l_sum = 0;
  if (p_f->m >= minLogSize && p_f->m <= maxLogSize)
  {
    l_numFactors = factorOffset(p_f->m + 1) - factorOffset(p_f->m);
    for (i = 0; i < l_numFactors; i++)
    {
      l_sum += fft_factors[factorOffset(p_f->m) + i];
    }
  }
  else
  {
    l_sum = p_f->n - 1;
  }
  return l_sum;
}

uint32_t numFactors(uint32_t m)
{
  uint32_t res;
  //One checks if one has precalculated factors for the desired fft size
  if (m >= minLogSize && m <= maxLogSize)
  {
    res = factorOffset(m + 1) - factorOffset(m);
    //l_ptr = fft_factors + factorOffset(m);
  }
  else res = 1;
  return res;
}

uint32_t* factors(uint32_t m)
{
  uint32_t * l_ptr = 0;
  //One checks if one has precalculated factors for the desired fft size
  if (m >= minLogSize && m <= maxLogSize)
  {
    l_ptr = fft_factors + factorOffset(m);
  }
  return l_ptr;
}

void DFT_direct(TfiniteField* p_f, uint32_t* pio_data, uint32_t* p_buffer)
{
  uint32_t i, j;
  uint32_t* log = p_f->log;
  uint32_t* exp = p_f->exp;
  uint32_t order = p_f->n - 1;
  for (i = 0; i < order; i++)
  {
    // p_buffer[i] will hold P(x**i)
    p_buffer[i] = 0;
    for (j = 0; j < order; j++)
    {
      uint32_t x = pio_data[j];
      if(x != 0)
      {
        p_buffer[i]^= exp[((uint32_t)log[x]+i*j)%order];
      }
    }
  }
  memcpy(pio_data, p_buffer, order * sizeof(uint32_t));
}

int multiplicative_FFT(TfiniteField* p_f, uint32_t* pio_data, uint32_t* p_buffer)
{
  uint32_t num_factors = numFactors(p_f->m);
  uint32_t* facts = factors(p_f->m);
  if(facts == NULL)
  {
    // we don't know the factoring of 2^m-1, hence we cannot proceed.
    return 1;
  }
  //multiplicative_FFT_with_factors(p_f, pio_data, p_buffer, facts, num_factors);
  const uint32_t* log = p_f->log;
  const uint32_t order = p_f->n - 1;
  for (uint32_t i = 0; i < order; i++) pio_data[i] = log[pio_data[i]];

  multiplicative_FFT_step_chaining(p_f,pio_data, p_buffer, facts, num_factors, order, 1);
  return 0;
}

void multiplicative_FFT_compute_step(
    TfiniteField* p_f,
    uint32_t* p_logdata,
    uint32_t* po_buffer,
    uint32_t p_local_exponent,
    uint32_t p_initial_exponent,
    uint32_t p_stride_glob)
{
  const uint32_t* exp = p_f->exp;
  const uint32_t logzero = p_f->logzero;
  const uint32_t order = p_f->n - 1;
  // note : parentheses below are REQUIRED to avoid overflows
  const uint32_t stride = p_stride_glob * (p_local_exponent / p_initial_exponent); // number of DFTs performed
  const uint32_t N = order / p_local_exponent; // size of each DFT
  const uint32_t buf_size = N * stride;
  {
    for (uint32_t i = 0; i < buf_size; i++) po_buffer[i] = 0;

    // perform 'stride' dfts of size N (direct computation)
    // input is in log form, result is in p_buffer and and in normal (not log) form
    uint32_t idx_i, idx_j;
    for (uint32_t i = 0; i < N; i++)
    {
      idx_j = 0;
      for (uint32_t j = 0; j < N; j++)
      {
        idx_i = stride * i;
        for (uint32_t s = 0; s < stride; s++)
        {
          //here,
          //idx_i = s + stride * i;
          //idx_j = s + stride * j;
          if (p_logdata[idx_j] != logzero)
          {
            uint32_t x = p_logdata[idx_j];
            po_buffer[idx_i] ^= exp[x];
            // at all times here, pio_data[idx_j] = d[idx_j] + p_local_exponent * i * j mod order
            // where d is the initial value of pio_data[idx_j]
            x += j * p_local_exponent;
            //x %= order;
            if (x >= order)  x -= order;
            p_logdata[idx_j] = x;
          }
          idx_i++;
          idx_j++;
        }
      }
    }
  }
}

void multiplicative_FFT_reorder_step(
    TfiniteField* p_f,
    uint32_t* p_buffer,
    uint32_t* po_logdata,
    uint32_t p_local_exponent,
    uint32_t p_initial_exponent,
    uint32_t p_stride_ext,
    uint32_t p_stride_glob)
{
  const uint32_t* log = p_f->log;
  const uint32_t logzero = p_f->logzero;
  const uint32_t order = p_f->n - 1;
  // note : parentheses below are REQUIRED to avoid overflows
  const uint32_t stride = p_stride_glob * (p_local_exponent / p_initial_exponent); // number of DFTs performed
  const uint32_t N = order / p_local_exponent; // size of each DFT
  const uint32_t group_size = p_stride_ext * p_stride_glob;
  const uint32_t m = p_stride_ext * p_initial_exponent;
  uint32_t idx_in = 0, idx_out = 0;
  for (uint32_t i = 0; i < N; i++)
  {
    uint32_t s_ext = 0, s_int = 0;
    idx_out = group_size * i;
    for (uint32_t s = 0; s < stride; s++)
    {
      //here,
      // s_ext = s % group_size
      // s_int = (s - s_ext) / group_size
      // idx_in = i * stride + s with s < stride
      // idx_out = s_ext + group_size * (i + s_int * N)
      // where s_ext < group_size and i < N
      // data is read at index idx_in, multiplied, then written back in output array
      // at position idx_out
      uint32_t x = log[p_buffer[idx_in]];
      if(x != logzero)
      {
        x += m * s_int * i;
        // s_int <= s / group_size < stride / group_size = p_local_exponent / p_initial_exponent / p_stride_ext
        // therefore m * s_int = p_stride_ext * p_initial_exponent * s_int < p_local_exponent
        // and i < N
        // hence p_stride_ext * p_initial_exponent * s_int * i < p_local_exponent * N = order
        // initially x < order
        // hence x %= order is equivalent to:
        if(x >= order) x -= order;
      }
      po_logdata[idx_out] = x;
      s_ext++;
      idx_in++;
      idx_out++;
      if(s_ext == group_size)
      {
        s_ext = 0;
        s_int++;
        idx_out = s_ext + group_size * (i + s_int * N);
      }
    }
  }
}

void multiplicative_FFT_step_chaining(
    TfiniteField* p_f,
    uint32_t* pio_data,
    uint32_t* p_buffer,
    const uint32_t* p_factors,
    uint32_t p_num_factors,
    uint32_t p_order,
    uint32_t e)
{
  uint32_t p_stride_ext = 1;
  for (uint32_t i = 0; i < p_num_factors; i++)
  {
    const uint32_t p_local_exponent = p_order / p_factors[i];
    multiplicative_FFT_compute_step(p_f, pio_data, p_buffer, p_local_exponent, e, e);
    if(i != p_num_factors - 1)
    {
      multiplicative_FFT_reorder_step(p_f, p_buffer, pio_data, p_local_exponent, e, p_stride_ext, e);
    }
    else
    {
      for (uint32_t j = 0; j < p_f->n-1; j++) pio_data[j] = p_buffer[j];
    }
    p_stride_ext *= p_factors[i];
  }
}

int evaluate_polynomial_multiplicative_FFT(
    TfiniteField* p_f,
    uint32_t* p_polynomial,
    uint32_t* p_eval,
    uint32_t p_numCoeffs,
    uint32_t* p_buffer)
{
  // if p_polynomial == p_eval :
  // prepare coeffs in p_buffer
  // use p_polynomial as buffer in chaining

  uint32_t* prepare_buf = p_polynomial == p_eval ? p_buffer : p_eval;
  uint32_t* aux_buf = p_polynomial == p_eval ? p_polynomial : p_buffer;

  uint32_t e = 1;
  const uint32_t order = p_f->n - 1;
  const uint32_t logzero = p_f->logzero;
  const uint32_t* log = p_f->log;
  const uint32_t* l_factors = factors(p_f->m);
  if(p_numCoeffs <= 1)
  {
    uint32_t value = p_numCoeffs == 1 ? p_polynomial[0] : 0;
    for(uint32_t i = 0; i < order; i++) p_eval[i] = value;
    return 0;
  }
  if(l_factors == 0)
  {
    // we don't know the factoring of 2^m-1, hence we cannot proceed.
    return 1;
  }

  uint32_t num_FFT_factors = numFactors(p_f->m);
  while (num_FFT_factors > 0 && e * l_factors[num_FFT_factors - 1] * p_numCoeffs <= order)
  {
    num_FFT_factors--;
    e *= l_factors[num_FFT_factors];
  }
  // the factors of e are the ones for which we won't have to compute a FFT.
  // this is why we start from the end in loop above (factors are sorted in ascending order).

  // compute the DFTs of polynomials Q_j(X) = sum_i (p_i x**ij) X**i, of order order/e, j = 0 ... e - 1
  // i.e. with base element x**e
  // j-th DFT outputs P(x**(j+e*k)), k = 0 ... order/e - 1.
  // piecing together these results for j = 0 ... e - 1
  // yields the evaluation of P on all x**k for k = 0 ... order - 1.

  // prepare all polynomials Q_i, in log form
  for(uint32_t i = 0; i < order / e; i++)
  {
    uint32_t x = i < p_numCoeffs ? log[p_polynomial[i]] : logzero;
    if(x == logzero) for(uint32_t j = 0; j < e; j++) prepare_buf[e * i + j] = logzero;
    else
    {
      for(uint32_t j = 0; j < e; j++)
      {
        // here, x = log(p_i) + i*j
        prepare_buf[e * i + j] = x;
        x += i;
        if(x >= order) x-= order;
      }
    }
  }
  // compute DFTs
  multiplicative_FFT_step_chaining(p_f, prepare_buf, aux_buf, l_factors, num_FFT_factors, order, e);
  return 0;
}
