#include "gf2_polynomials_misc.h"

#define g_wordSize (8 * sizeof(unsigned long))

uint32_t word_size()
{
  return g_wordSize;
}

size_t wordsize_to_max_degree(uint32_t p_word_size)
{
  return p_word_size * g_wordSize - 1;
}

size_t degree_packed(unsigned long* p_polynomial, size_t p_max_degree)
{
  size_t wordIdx, offset;
  uint32_t maxOffset;
  if(p_max_degree == -1u) return -1u;
  degree_to_packed_pos(p_max_degree, &wordIdx, &maxOffset);

  while(p_polynomial[wordIdx] == 0 && wordIdx > 0)
  {
    wordIdx--;
  }
  if(p_polynomial[wordIdx] == 0)
  {
    // necessarily wordIdx == 0, and polynomial is 0
    return -1u;
  }
  offset = g_wordSize - 1;
  while((p_polynomial[wordIdx] >> offset) == 0)
  {
    offset--;
  }
  return g_wordSize * wordIdx + offset;
}

uint32_t degree_packed_u32(unsigned long* p_polynomial, uint32_t p_max_degree)
{
  return (uint32_t) degree_packed(p_polynomial, (size_t) p_max_degree);
}

size_t degree_to_packed_size(size_t p_degree)
{
  if(p_degree == -1u) return 0;
  return (p_degree / g_wordSize) + 1;
}

uint32_t degree_to_packed_size_u32(uint32_t p_degree)
{
  if(p_degree == -1u) return 0;
  return (p_degree / g_wordSize) + 1;
}

inline void degree_to_packed_pos(size_t p_degree, size_t* po_word_idx, uint32_t* po_word_offset)
{
  *po_word_idx = p_degree / g_wordSize;
  if(po_word_offset) *po_word_offset = (uint32_t) (p_degree - g_wordSize * *po_word_idx);
}

inline void degree_to_packed_pos_u32(uint32_t  p_degree, uint32_t* po_word_idx, uint32_t* po_word_offset)
{
  *po_word_idx = p_degree / g_wordSize;
  if(po_word_offset) *po_word_offset = (uint32_t) (p_degree - g_wordSize * *po_word_idx);
}

void set_bit(unsigned long* p_polynomial, size_t p_degree, uint32_t p_val)
{
  size_t word_idx;
  uint32_t word_offset;
  degree_to_packed_pos(p_degree, &word_idx, &word_offset);
  if(p_val)
  {
    p_polynomial[word_idx] |= (1ul << word_offset);
  }
  else
  {
    p_polynomial[word_idx] &= ~(1ul << word_offset);
  }
}

uint32_t get_bit(unsigned long* p_polynomial, size_t p_degree)
{
  size_t wordIdx;
  uint32_t wordOffset, res;
  degree_to_packed_pos(p_degree, &wordIdx, &wordOffset);
  res = (p_polynomial[wordIdx] >> wordOffset) & 1;
  return res;
}

void add_bit(unsigned long* p_polynomial, size_t p_degree)
{
  size_t wordIdx;
  uint32_t wordOffset;
  degree_to_packed_pos(p_degree, &wordIdx, &wordOffset);
  p_polynomial[wordIdx] ^= (1ul << wordOffset);
}




