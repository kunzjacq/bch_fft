#pragma once

#include <stdint.h>

#include "finite_field.h"

#ifdef	__cplusplus
extern "C" {
#endif

uint32_t numFactors(uint32_t m);
uint32_t* factors(uint32_t m);
uint32_t factorOffset(uint32_t p_log);
uint32_t getFactorSum(TfiniteField* p_f);


/**
 * @brief DFT_direct
 * Fourier transform of input polynomial P defined by 'pio_data', i.e. evaluation of P(x)
 * for all x in finite field f defined by 'p_f'. Output is in-place.
 * @param p_f
 * @param pio_data
 * @param p_buffer
 * a buffer of (p_f->n - 1) * 4 bytes.
 */

void DFT_direct(
    TfiniteField* p_f,
    uint32_t* pio_data,
    uint32_t* p_buffer);

/**
 * @brief multiplicative_FFT
 * Performs a FFT of an input polynomial P described by 'pio_data' with coefficients in finite
 * field 'p_f', using canonical primitive element x of p_f, i.e. outputs
 * P(x^i), i = 0 ... order - 1, order = p_f->n - 1.
 * Uses the factorization of order to speed up the computation if this factorization is present
 * in a precomputed table.
 * @param p_f
 * The binary finite field description.
 * @param pio_data
 * The array of size 'order' for input polynomial and output DFT result.
 * @param p_buffer
 * a buffer of size 'order'.
 */

int multiplicative_FFT(TfiniteField* p_f, uint32_t* pio_data, uint32_t* p_buffer);

/**
 * @brief multiplicative_FFT_compute_step
 * computes stride = 'p_stride_glob' * 'p_local_exponent' / 'p_initial_exponent' DFTs in
 * binary field 'p_f'.
 * DFTs are of order N = order / 'p_local_exponent', order = (p_f->n - 1),
 * and are computed using direct formula (with Chien trick), on data D = 'pio_data'.
 * The 'stride' groups of data are interleaved:
 * group 0 = [D[0], D[stride], ..., D[stride * (N-1)]],
 * group 1 = [D[1], D[stride+1], ..., D[stride * (N-1) + 1]],
 * ...
 * group 'stride - 1' = [D[stride - 1], D[2* stride-1], ..., D[stride * N - 1]];
 * the FFT is computed for element y = x^'p_local_exponent' where x is the primitive element whose
 * representation is 2.
 * buffer size = stride * N = p_stride_glob * order / p_initial_exponent.
 *
 * It is assumed that p_initial_exponent | p_local_exponent and
 *  p_local_exponent | order.
 * @param p_f
 * the finite field description.
 * @param p_logdata
 * the input data, in log form.
 * @param po_buffer
 * the output data, in plain form.
 * @param p_local_exponent
 * see general description.
 * @param p_initial_exponent
 * see general description.
 * @param p_stride_glob
 * see general description.
 */

void multiplicative_FFT_compute_step(
    TfiniteField *p_f,
    uint32_t* p_logdata,
    uint32_t* po_buffer,
    uint32_t p_exponent,
    uint32_t p_initial_exponent,
    uint32_t p_stride_glob);


/**
 * @brief multiplicative_FFT_reorder_step
 * shuffles and multiplies elements on input, which is composed of 'stride' interleaved series of
 * size N = order / 'p_local_exponent', order = (p_f->n - 1) as with multiplicative_FFT_compute_step.
 * with group_size = p_stride_ext * p_stride_glob, the 'stride' series are regrouped in
 * stride / group_size subgroups of size group_size.
 * Within a subgroup, series are interleaved as on input.
 * Within subgroup if index s_int, each i-th term is multiplied by
 * i * s_int * p_stride_ext * p_initial_exponent.
 *
 * input is in log form, output in plain form.
 * @param p_f
 * the finite field used.
 * @param p_buffer
 * input data, in plain form.
 * @param po_logdata
 * output data, in log form.
 * @param p_local_exponent
 * see general description.
 * @param p_initial_exponent
 * see general description.
 * @param p_stride_ext
 * a substride for reordering of data.
 * it is assumed that p_stride_ext * p_stride_glob | stride = p_local_exponent / p_initial_exponent
 * @param p_stride_glob
 * a substride for reordering of data.
 * it is assumed that p_stride_ext * p_stride_glob | stride = p_local_exponent / p_initial_exponent
 * p_stride_glob plays no role in post-multiplication, whereas p_stride_ext does.
 */

void multiplicative_FFT_reorder_step(
    TfiniteField* p_f,
    uint32_t* p_buffer,
    uint32_t* po_logdata,
    uint32_t p_local_exponent,
    uint32_t p_initial_exponent,
    uint32_t p_stride_ext,
    uint32_t p_stride_glob);

/**
 * @brief multiplicative_FFT_step_chaining
 * correctly chain calls to multiplicative_FFT_compute_step and multiplicative_FFT_reorder_step
 * to compute one or several interleaved fourier transforms of an order with known factors.
 *
 * to compute p_stride_glob DFT for element x^'initial_exponent', of order 'order' / 'initial_exponent':
 * with factors[0] * ... * factors[k - 1] = order, where factors[i] is prime
 * (and initial_exponent | order, initial_exponent < order)
 * perform k calls to multiplicative_FFT_compute_step, with i = 0 ... k-1,
 * with a call to multiplicative_FFT_reorder_step between these calls.
 * local_exponent = order / p_factors[i] so that initial_exponent | local_exponent
 * stride = p_stride_glob * local_exponent / initial_exponent
 * stride_ext = factors[0] * ... * factors[i-1]
 * N = p_factors[i]
 *
 * to understand the algorithm, it is enough to understand the case 'initial_exponent' = 1
 * since the general case is the same as this case with x^'initial_exponent' replacing x;
 * and p_stride_glob = 1, since p_stride_glob > 1 performs p_stride_glob executions of the
 * algorithm with p_stride_glob = 1 on interleaved input and output.
 * @param p_f
 * @param pio_data
 * @param p_buffer
 * @param p_factors
 * @param p_num_factors
 * @param p_order
 * @param e
 */

void multiplicative_FFT_step_chaining(
    TfiniteField* p_f,
    uint32_t* pio_data,
    uint32_t* p_buffer,
    const uint32_t *p_factors,
    uint32_t p_num_factors,
    uint32_t order,
    uint32_t e);

/**
 * @brief evaluate_poynomial_multiplicative_FFT
 * evaluates the input polynomial P on all field values, using DFT. takes into account the
 * polynomial degree and the factorization of the finite field multiplicative group order
 * to optimize the computation.
 * output is in p_result: p_result[i] = P(x^i) where x is the canonical primitive element of
 * the finite field.
 * @param p_f
 * finite field description
 * @param p_polynomial
 * array of polynomial coefficients. The input polynomial is not modified during computation.
 * @param p_eval
 * result buffer. buffer can be equal to p_polynomial, or disjoint from it, but a non-trivial
 * overlap is forbidden and will lead to incorrect results at best.
 * @param p_numCoeffs
 * number of polynomial coefficients (i.e. degree + 1)
 * @param p_buffer
 * a buffer of size p_f->n - 1 (disjoint from all other buffers).
 */

int evaluate_polynomial_multiplicative_FFT(
    TfiniteField* p_f,
    uint32_t* p_polynomial,
    uint32_t* p_eval,
    uint32_t p_numCoeffs,
    uint32_t *p_buffer);

#ifdef	__cplusplus
}
#endif


