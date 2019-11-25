#include <stdio.h>
#include <string.h> // for memset / memcpy
#include <assert.h>

#include "gf2_polynomials_misc.h"
#include "gf2_polynomials.h"

static uint32_t g_wordSize = (8 * sizeof (unsigned long));

void euclidean_division_packed_inplace_gf2(
    unsigned long* pio_reductee,
    uint32_t p_reductee_degree,
    unsigned long* p_reductor,
    uint32_t p_reductor_degree,
    unsigned long* po_dividend)
{
  uint32_t i, uPrime, uFinal, offset, u;
  uint32_t reductorNumWords;
  unsigned long bit, word1, word2;

  reductorNumWords = degree_to_packed_size_u32(p_reductor_degree);
  if (po_dividend)
  {
    memset(po_dividend, 0, degree_to_packed_size(p_reductee_degree - p_reductor_degree) * sizeof (unsigned long));
  }
  if(p_reductor_degree > p_reductee_degree) return;
  for (i = p_reductee_degree; i != p_reductor_degree - 1; i--)
  {
    // introduce a new monomial of the reductee
    bit = get_bit(pio_reductee, i);
    if (bit)
    {
      // replace the term X^i by its reduction
      // by adding [reductor times X^(i-p_reductorDegreed)] to the reductee
      degree_to_packed_pos_u32(i - p_reductor_degree, &u, &offset);
      if (po_dividend)
      {
        po_dividend[u] ^= 1ul << offset;
      }
      if (offset == 0)
      {
        for (uPrime = 0; uPrime < reductorNumWords; uPrime++)
        {
          pio_reductee[u + uPrime] ^= p_reductor[uPrime];
        }
      }
      else
      {
        for (uPrime = 0; uPrime < reductorNumWords - 1; uPrime++)
        {
          word1 = p_reductor[uPrime] << offset;
          word2 = p_reductor[uPrime] >> (g_wordSize - offset);
          uFinal = u + uPrime;
          pio_reductee[uFinal] ^= word1;
          pio_reductee[uFinal + 1] ^= word2;
        }
        // the last word is treated as a special case to avoid overflows
        uPrime = reductorNumWords - 1;
        uFinal = u + uPrime;
        word1 = p_reductor[uPrime] << offset;
        pio_reductee[uFinal] ^= word1;
        if (degree_to_packed_size(p_reductee_degree) > uFinal + 1)
        {
          word2 = p_reductor[uPrime] >> (g_wordSize - offset);
          pio_reductee[uFinal + 1] ^= word2;
        }
      }
    }
  }
}
