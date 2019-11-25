#pragma once

#include <stdint.h>
#include <malloc.h>
#include <string.h>

#ifdef	__cplusplus
extern "C" {
#endif

#define g_maxDegree 30

typedef struct
{
  // log of field size
  uint32_t m;
  // field size
  uint32_t n;
  // a conventional value s.t.
  // log(0) = logzero and exp[log[0]] = 0
  uint32_t logzero;
  // array i -> x**i
  uint32_t* exp;
  // array x**i -> i
  uint32_t* log;
} TfiniteField;

/**
 * @brief create_binary_finite_field
 * creates a structure describing a finite field of characteristic 2, of size 2**p_degree.
 * the internal representation of the field is a logarithm and exponent table for a primitive
 * element whose representation is 2. multiplication is performed with these tables, addition
 * is the xor operation.
 * This approach only works for small finite fields. here it is made available for sizes from
 * 2**2 to 2**30, inclusive.
 * @param p_degree
 * the finite field degree over F_2.
 * @return
 * a TfiniteField structure.
 */
TfiniteField create_binary_finite_field(uint32_t p_degree);

/**
 * @brief free_finite_field
 * frees a data structure describing a finite field.
 * @param p_f
 * the finite field to process.
 */

void free_finite_field(TfiniteField * p_f);


/**
 * @brief inverse
 * returns 1/a if a is nozero, 0 otherwise
 * @param p_f
 * the finite field used.
 * @param p_a
 * value of a
 * @return
 * 1 / a
 */

uint32_t inverse(TfiniteField * p_f, uint32_t p_a);


/**
 * @brief square
 * squares a.
 * @param p_f
 * the finite field used.
 * @param p_a
 * value of a
 * @return
 * a * a
 */

uint32_t square(TfiniteField * p_f, uint32_t p_a);

/**
 * @brief multiply
 * multiplies a by b.
 * @param p_f
 * the finite field used.
 * @param p_a
 * value of a
 * @param p_b
 * value of b
 * @return
 * a * b
 */

uint32_t multiply(TfiniteField * p_f, uint32_t p_a, uint32_t p_b);

/**
  @brief multiply_by_log
 * returns a * b, where a is given by p_a and b is given by its log p_bLog.
 * b is assumed to be non zero (i.e. p_bLog != p_f->logzero)
 * @param p_f
 * the finite field used
 * @param p_a
 * value of a
 * @param p_b_log
 * value of b (nonzero), in log form.
 * @return
 * a * b
 */

uint32_t multiply_by_log(TfiniteField * p_f, uint32_t p_a, uint32_t p_b_log);

#ifdef	__cplusplus
}
#endif


