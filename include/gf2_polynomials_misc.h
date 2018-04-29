#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C"
{
#endif

// the typedef for non-packed values in GF(2)
typedef unsigned char binary_coeff;

// dense packed values are unsigned longs, that hold 8 * sizeof(unsigned long) coefficients

// returns the number of coefficients stored in a word (unsigned long)
uint32_t word_size();

// returns the maximum degree of a polynomial stored in 'p_word_size' words
size_t wordsize_to_max_degree(uint32_t p_word_size);

// degree of a packed polynomial 'p_polynomial' of maximum degree 'p_max_degree'
size_t degree_packed(unsigned long* p_polynomial, size_t p_max_degree);

// same as above, with uint32_t
uint32_t degree_packed_u32(unsigned long* p_polynomial, uint32_t p_max_degree);

// size in number of words of a packed polynomial of degree 'p_degree'
size_t degree_to_packed_size(size_t p_degree);

// same as above, with uint32_t
uint32_t degree_to_packed_size_u32(uint32_t p_degree);

// compute the position of a monomial of a given degree in a packed polynomial
void degree_to_packed_pos(size_t p_degree, size_t* po_word_idx, uint32_t* po_word_offset);

void degree_to_packed_pos_u32(uint32_t  p_degree, uint32_t* po_word_idx, uint32_t* po_word_offset);

// set to 'p_val' coefficient of monomial of degree 'p_degree' in polynomial 'p_polynomial'
void set_bit(unsigned long* p_polynomial, size_t p_degree, uint32_t p_val);

// get coefficient of monomial of degree 'p_degree' in polynomial 'p_polynomial'
uint32_t get_bit(unsigned long* p_polynomial, size_t p_degree);

// add 1 to coefficient of monomial of degree 'p_degree' in polynomial 'p_polynomial'
void add_bit(unsigned long* p_polynomial, size_t p_degree);


#ifdef __cplusplus
}
#endif
