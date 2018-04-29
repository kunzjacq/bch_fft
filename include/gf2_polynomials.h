#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @brief euclidean_division_packed_inplace_gf2
 * if P and Q polynomials with coefficients of GF(2), euclidean_division_packed_inplace_gf2
 * computes polynomials R, and if desired, S, s.t. P = Q S + R and (R = 0 or deg(R) < deg(Q)),
 * that is, computes the Euclidean division of P by Q.
 * P is replaced by R in the process of computation, that is, R is computed in- place.
 * All polynomials are dense and packed:
 * coefficient of degree i of a polynomial stored in an unsigned long array p is the
 * (i mod k) - th bit of p[i/k],  if unsigned longs are k-bit values.
 * @param pio_reductee
 * the polynomial P, that contains R on exit
 * @param p_reductee_degree
 * degree of P
 * @param p_reductor
 * the polynomial Q
 * @param p_reductor_degree
 * degree of Q
 * @param po_dividend
 * the address of an array where S is output, of size sufficient to hold a polynomial of
 * degree deg(P) - deg(Q)
 * (on output, only the coefficients up to that degree are guaranteed to be correct.)
 * if po_dividend is 0, S is not computed.
 */

void euclidean_division_packed_inplace_gf2(
    unsigned long* pio_reductee,
    uint32_t p_reductee_degree,
    unsigned long* p_reductor,
    uint32_t p_reductor_degree,
    unsigned long* po_dividend);

#ifdef __cplusplus
}
#endif
