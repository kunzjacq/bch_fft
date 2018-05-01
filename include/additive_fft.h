#pragma once

#include "finite_field.h"
#include "gf2_extension_polynomials.h"

// wrapper around one of the two functions below
// buffers must be distinct ; source buffer is modified
void evaluate_polynomial_additive_FFT(
    TfiniteField* p_f,
    uint32_t* p_poly,
    uint32_t* po_result,
    uint32_t p_num_terms);

/*
 * additive_DFT, additive_DFT_opt
 * in-place Von zur Gathen-Gerhard additive FFT of input polynomial
 * explained for instance in Todd Mateer PhD thesis
 * available at
 * https://tigerprints.clemson.edu/all_dissertations/231/
 * or
 * https://cr.yp.to/f2mult.html
 *
 * outputs p(i) for all i in a binary finite field.
 * if field size is 2^n, performs at O(n (log n)^2) additions and multiplications in the field.
 *
 * note: the same PhD thesis introduces a faster algorithm (chap 3.5); see the second link above.
 */

/**
 * @brief additive_DFT_opt
 * @param p_f
 * finite field used
 * @param p_poly
 * input polynomial to evaluate and output result
 * @param p_poly_degree
 * degree of input polynomial p_poly, or 0. If > 0, unnecessary steps of the
 * algorithm are skipped.
 */

void additive_DFT_opt(TfiniteField* p_f, uint32_t* p_poly, uint32_t p_poly_degree);

/**
 * @brief additive_DFT
 * @param p_f
 * finite field used
 * @param p_poly
 * input polynomial to evaluate and output result
 * @param p_poly_degree
 * degree of input polynomial p_poly, or 0. If > 0, unnecessary steps of the
 * algorithm are skipped.
 */

void additive_DFT(TfiniteField* p_f, uint32_t* p_poly, uint32_t p_poly_degree);
