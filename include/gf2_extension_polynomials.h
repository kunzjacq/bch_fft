#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

#include "finite_field.h"

/**
 * @brief evaluate_sparse_polynomial
 * evaluates a sparse polynomial at one point.
 *
 * @param p_f
 * finite field definition
 * @param p_coeffs
 * monomial coefficients.
 * @param p_degrees
 * monomial degrees.
 * @param p_num_coeffs
 * number of monomials.
 * @param p_y
 * point where the polynomial is evaluated.
 * @return
 * the value of the polynomial.
 */

uint32_t evaluate_sparse_polynomial(
    TfiniteField* p_f,
    uint32_t* p_coeffs,
    uint32_t* p_degrees,
    uint32_t p_num_monomials,
    uint32_t p_y);

/**
 * @brief evaluate_sparse_polynomial_log
 * evaluates a sparse polynomial at one point; coefficients of the polynomial are given by their logs.
 * @param p_f
 * finite field definition
 * @param p_log_coeffs
 * logs of monomial coefficients.
 * @param p_degrees
 * degrees corresponding to coefficients of the polynomial.
 * @param p_num_coeffs
 * total number of monomials.
 * @param p_y
 * point where the polynomial is evaluated.
 * @return
 * the value of the polynomial.
 */

uint32_t evaluate_sparse_polynomial_log(
    TfiniteField* p_f,
    uint32_t* p_log_coeffs,
    uint32_t* p_degrees,
    uint32_t num_coeffs,
    uint32_t p_y);


/**
 * @brief multiply_polynomial_by_monomial
 * multiply in-place some polynomial 'pio_poly' of degree <= 'p_degree' by
 * (X - x**'p_nonzeroRootLog')
 * @param p_f
 * the finite field used, of primitive element x with representation 2
 * @param pio_poly
 * the input and output polynomial. buffer must be able to hold at least p_degree + 2 coefficients.
 * @param p_degree
 * see general description
 * @param p_nonzero_root_log
 * see general description
 */

void multiply_polynomial_by_monomial(
    TfiniteField* p_f,
    uint32_t* pio_poly,
    uint32_t p_degree,
    uint32_t p_nonzero_root_log);

/**
 * @brief euclidean_division_sparse_reductor
 * euclidean division of a dense polynomial by a sparse polynomial in a binary field.
 * @param p_f
 * finite field description
 * @param p_reductee
 * coefficients of polynomial to reduce.
 * @param p_reductee_degree
 * degree of polynomial to reduce.
 * @param p_reductor_coeffs
 * reductor monomial coefficients
 * @param p_reductor_degrees
 * reductor monomial degrees
 * @param p_reductor_num_terms
 * number of terms in reductor.
 * @param po_remainder
 * output array for remainder, of size (reductor degree).
 * @param po_dividend
 * output array for dividend, or 0. The dividend is written if po_dividend != 0.
 * if this array is nonzero, it must be of size p_reductee_degree - (reductor degree) + 1.
 * @param p_dividend_factor_log
 * the dividend is multiplied by this value in po_dividend.
 */

void euclidean_division_sparse_reductor(
    TfiniteField* p_f,
    const uint32_t *p_reductee,
    uint32_t p_reductee_degree,
    uint32_t* p_reductor_coeffs,
    uint32_t* p_reductor_degrees,
    uint32_t p_reductor_num_terms,
    uint32_t* po_remainder,
    uint32_t* po_dividend,
    uint32_t p_dividend_factor_log);

#ifdef __cplusplus
}
#endif
