#pragma once

#include "finite_field.h"

#ifdef	__cplusplus
extern "C"
{
#endif
  
typedef struct _TbchCode
{
  TfiniteField m_field;
  // code minimum distance
  uint32_t m_distance;
  // the generating polynomial, in packed form
  unsigned long* m_packed_gen_poly;
  // degree of said polynomial
  uint32_t m_gen_poly_degree;
  // code word length
  uint32_t m_length;
} TbchCode;

// a type for a buffer for BCH decoding
typedef struct _TbchCodeBuffer
{
  // buffer of size p_codePtr->m_field.n
  uint32_t* m_syndrome;
  // buffer of size p_codePtr->m_field.n
  uint32_t* m_buffer;
  // buffers for polynomials used by berlekamp-massey algorithm.
  // m_A has size p_codePtr->m_field.n, to serve other buffering purposes outside of B-M.
  // other buffers have size p_codePtr->m_distance + 1.
  uint32_t* m_A;
  uint32_t* m_B;
  // pointer to the error location polynomial.
  uint32_t* m_elp;
  // buffer storing location of errors.
  uint32_t* m_error;
} TbchCodeBuffer;

/**
 * @brief bch_gen_code
 * Generate a BCH code of codeword size p_codeword_size and of minimum distance
 * p_minimum_distance
 * @param p_codeword_size
 * @param p_minimum_distance
 * the desired minimum distance for the generated code. (=2*c+1 where c is the desired
 * error-correcting capability).
 * @return
 * the generated code
 */

TbchCode bch_gen_code(uint32_t p_codeword_size, uint32_t p_minimum_distance);

/**
 * @brief bch_build_polynomial
 * Compute the generator polynomial of a binary BCH code. Fist generate the
 * cycle sets modulo 2**m - 1, cycle[][] =  (i, 2*i, 4*i, ..., 2^l*i). Then
 * determine those cycle sets that contain integers in the set of (d-1)
 * consecutive integers {1..(d-1)}. The generator polynomial is calculated
 * as the product of linear factors of the form (X+x**i), for every i in
 * the above cycle sets.
 * @param p_f
 * the finite field used. note: it is copied into the output structure, along with the
 * pointer it contains. The TbchCode structure then 'owns' the finite field and destroys it
 * when it is destroyed, hence the FF should not be destroyed itself after calling
 * generateCodePolynomial.
 * @param p_minimumDistance
 * code minimum distance.
 * @return
 * a TbchCode structure.
 */
TbchCode bch_build_polynomial(TfiniteField* p_f, uint32_t p_minimumDistance);

/**
 * @brief bch_explore_cycles
 * forms the set S' of all conjugates of elements
 * in S = {p_offset, p_offset + 1, ..., p_offset + p_numValues - 1}
 * modulo p_sz, with p_sz having form 2**m - 1.
 * the set has maximum size m * p_numValues, since each element is S has at most m conjugates.
 * if p_zerosPtr != 0, returns elements of S' sequentially in p_zerosPtr[0], [1],...
 * else returns only the size of S
 * @param p_f
 * the finite field used.
 * @param p_numValues: description of S, see above
 * @param p_offset: description of S, see above
 * @param p_zerosPtr 0 or a pointer to the result array, of size at least
 *   min(2^m - 1, m * p_minimumDistance)
 * @param p_cyclePtr a buffer of size m
 * @param p_remainingPtr a buffer of size 2**m = p_sz + 1
 * @param p_binaryPolynomials an array of unsigned longs of size p_numValues
 * @return the size of S
 */

void bch_explore_cycles(
    TfiniteField *p_f, uint32_t p_minimumDistance, uint32_t p_offset,
    uint32_t* p_remainingPtr, unsigned long *p_binaryPolynomials, uint32_t *po_numPolys,
    uint32_t *po_totalDegree);

/**
 * @brief bch_reduce_syndrome
 * reduces 'pio_encoded_data' by 'p_code->m_generating_polynomial'.
 * @param p_code
 * the error correcting code description.
 * @param pio_encodedData
 * input and output GF(2) polynomial, in packed form.
 */

void bch_reduce_syndrome(TbchCode *p_code, unsigned long* pio_encoded_data);

/**
 * @brief bch_eval_syndrome
 * Evaluate syndrome (as a polynomial) on points 0, ..., p_code->m_distance. outputs 1 if there is
 * at least one error, 0 otherwise.
 * Uses FFT or direct evaluation depending of p_code->m_distance. Pre-reduces syndrome by
 * generating polynomial if using direct evaluation.
 * @param p_code
 * the code description
 * @param p_buffer
 * the buffer structure. syndrome evaluation is stored in p_buffer-> m_syndrome.
 * @param p_packedData
 * the input syndrome, in packed form.
 * @return
 * 1 if there are errors to correct, 0 otherwise.
 */

int32_t bch_eval_syndrome(TbchCode* p_code, TbchCodeBuffer *p_buffer, unsigned long* p_packed_syndrome);

/**
 * @brief bch_compute_elp
 * compute the error location polynomial from syndrome with Berlekamp-Massey algorithm.
 * syndrome is in p_buffer->m_syndrome, B-M buffers are in p_buffer->m_{A,B,C}.
 * @param p_codePtr
 * @param p_buffer
 * @return
 * the degree of the ELP. The ELP itself is in p_buffer->m_elp
 * (which points to one of the buffers p_buffer->m_{A,B,C})
 */

uint32_t bch_compute_elp(TbchCode* p_code, TbchCodeBuffer* p_buffer);

/**
 * @brief bch_locate_errors
 * locate errors which are roots of the ELP, by evaluating the ELP using either a FFT or a direct
 * computation. Returns the number of roots found.
 * @param p_code
 * the error-correcting code description.
 * @param p_buffer
 * the buffer structure. ELP is in p_buffer->m_elp. auxiliary buffers are used.
 * the roots found are listed in p_buffer-> m_errors.
 * @param p_elp_degree
 * the ELP degree.
 * @return
 * the number of roots found.
 */

uint32_t bch_locate_errors(TbchCode* p_code, TbchCodeBuffer *p_buffer, uint32_t p_elp_degree);

/**
 * @brief bch_decode
 * chain all decoding operations.
 * @param p_code
 * the bch code description.
 * @param pio_packedData
 * the input codeword, and output corrected codeword.
 * @param p_buffer
 * the bch buffer structure.
 * @param po_corrected
 * if this pointer is non-zero, the number of corrections made is returned in this variable.
 * @return
 */

uint32_t bch_decode(TbchCode* p_code, unsigned long* pio_packed_data,
    TbchCodeBuffer* p_buffer, uint32_t *po_corrected);

/**
 * @brief bch_encode
 * Encode a block of data into a valid BCH codeword.
 * @param p_code
 * the BCH code description.
 * @param po_codeword
 * the output codeword.
 * @param p_data
 * the input data.
 * @param p_data_size
 * the input data size.
 * @return
 * 0 on success (if the data blocks fit into a codeword), 1 on error.
 */

int bch_encode(TbchCode *p_code, unsigned long* po_codeword, uint32_t* p_data, uint32_t p_data_size);

/**
  @brief bch_find_best_offset
 * Finds the offset that minimizes the degree of the generating polynomial, in a given
 * finite field p_f and for a distance p_minimum_distance.
 * For offset c,the generating polynomial has roots
 * x**c, ..., x**(c + p_minimum_distance-1), and their conjugates
 * where x is the canonical primitive element of the finite field.
 * In practice one never finds a better offset than 0, hence this function is not used;
 * all other functions assume that the generating polynomial has roots
 * x**0, ..., x**(p_minimum_distance-1).
 * @param p_f
 * the finite field used.
 * @param p_minimum_distance
 * See general description.
 */

void bch_find_best_offset(TfiniteField* p_f, uint32_t p_minimum_distance);

// misc memory management functions
TbchCodeBuffer bch_alloc_buffers(TbchCode* p_code);
void bch_free_code(TbchCode *p_code);
void bch_free_buffers(TbchCodeBuffer* p_buffer);

#ifdef	__cplusplus
}
#endif
