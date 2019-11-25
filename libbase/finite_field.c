#include <stdlib.h>

#include "finite_field.h"

TfiniteField create_binary_finite_field(uint32_t p_degree)
{
  static const uint32_t min_degree = 2;
  static const uint32_t max_degree = g_maxDegree;
  static const uint32_t monomial_degree[] =
  {
    1, 1, 1, 2, 1, 1, 4, 5, 6, 4, 3, 2, 3, 4, 7, 1, 3, 4, 1,
    11, 12, 1, 2, 3, 5, 3, 7, 1, 5, 6, 3, 2, 1, 5, 1, 3, 4,
    3, 1, 3, 4, 1, 2, 5, 1, 2, 1
  };
  static uint32_t offset[] = {
    0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 15, 18, 21, 22, 25, 26,
    27, 30, 31, 32, 33, 34, 37, 38, 41, 42, 43, 44, 45, 46, 47
  };
  TfiniteField f;
  f.log = NULL;
  f.exp = NULL;
  f.m = p_degree;
  f.n = (1 << p_degree);
  f.logzero = f.n - 1;

  if (f.m < min_degree || f.m > max_degree)
  return f; // we don't have a polynomial for the requested size

  f.exp = (uint32_t*) calloc(f.n, sizeof (uint32_t));
  f.log = (uint32_t*) calloc(f.n, sizeof (uint32_t));
  f.exp[f.logzero] = 0;
  const uint32_t d = f.m - min_degree;
  uint32_t poly = 1;
  for (uint32_t i = offset[d]; i < offset[d + 1]; i++) poly ^= 1u << monomial_degree[i];
  poly ^= f.n;
  uint32_t c = 1;
  for (uint32_t i = 0; i < f.n - 1; i++)
  {
    // here c = x**i where i is the chosen primitive element of the finite field
    // whose representation is 2
    f.exp[i] = c;
    c <<= 1;
    if (c & f.n) c ^= poly;
  }
  // build log table
  f.log[0] = f.logzero;
  for (uint32_t i = 0; i < f.n - 1; i++) f.log[f.exp[i]] = i;
  return f;
}

void free_finite_field(TfiniteField * p_f)
{
  free(p_f->exp);
  p_f->exp = 0;
  free(p_f->log);
  p_f->log = 0;
}

uint32_t inverse(TfiniteField * p_f, uint32_t p_a)
{
  uint32_t order = p_f->n - 1;
  if (p_a == 0) return 0;
  if (p_a == 1) return 1;
  uint32_t loginv = order - p_f->log[p_a];
  return p_f->exp[loginv];
}

uint32_t square(TfiniteField * p_f, uint32_t p_a)
{
  uint32_t order = p_f->n - 1;
  if (p_a == 0)  return 0;
  uint32_t logsq = 2 * p_f->log[p_a];
  if (logsq >= order) logsq -= order;
  return p_f->exp[logsq];
}

uint32_t multiply(TfiniteField * p_f, uint32_t p_a, uint32_t p_b)
{
  uint32_t order = p_f->n - 1;
  if (p_a == 0 || p_b == 0)  return 0;
  uint32_t logsum = p_f->log[p_a] + p_f->log[p_b];
  if (logsum >= order) logsum -= order;
  return p_f->exp[logsum];
}

uint32_t multiply_by_log(TfiniteField * p_f, uint32_t p_a, uint32_t p_b_log)
{
  uint32_t order = p_f->n - 1;
  if (p_a == 0) return 0;
  uint32_t logsum = p_f->log[p_a] + p_b_log;
  if (logsum >= order) logsum -= order;
  return p_f->exp[logsum];
}
