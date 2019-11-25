#include <stdio.h>
#include <string.h> // for memset etc

#include "gf2_extension_polynomials.h"

uint32_t evaluate_sparse_polynomial(
    TfiniteField* p_f,
    uint32_t* p_coeffs,
    uint32_t* p_degrees,
    uint32_t p_num_monomials,
    uint32_t p_y)
{
  uint32_t j,res;
  uint32_t order =  order = p_f->n-1;
  if(p_y == 0 && p_num_monomials > 0) res = p_degrees[0] == 0 ? p_coeffs[0] : 0;
  else
  {
    res = 0;
    uint64_t v;
    for(j = 0; j < p_num_monomials; j++)
    {
      if(p_coeffs[j])
      {
        // we want to compute (p_fieldPtr->log[i] * p_polyDegrees[j]) % (p_fieldPtr->n-1)
        // but a direct computation with uint32_t may overfow
        v = (uint64_t) p_f->log[p_y];
        v *= p_degrees[j];
        v = v % order;
        uint32_t s = multiply_by_log(p_f, p_coeffs[j], (uint32_t) v);
        res^= s;
      }
    }
  }
  return res;
}

static uint32_t power(
    TfiniteField* p_f,
    uint32_t p_val,
    uint32_t p_exp)
{
  if(p_exp == 0) return 1;
  int32_t pos = 31;
  uint32_t res = p_val;
  while((pos > 0) && (((p_exp >> pos) & 1) == 0)) pos--;
  for(; pos != 0; pos--)
  {
    res = square(p_f, res);
    if((p_exp >> (pos-1)) & 1) res = multiply(p_f, res, p_val);
  }
  return res;
}

uint32_t evaluate_sparse_polynomial_ref(
    TfiniteField* p_f,
    uint32_t* p_coeffs,
    uint32_t* p_degrees,
    uint32_t p_num_monomials,
    uint32_t p_y)
{
  uint32_t j, res;
  uint32_t curr_exp = 0;
  uint32_t curr_power = 1; // = p_y**curr_exp at all times
  if(p_y == 0 && p_num_monomials > 0) res = p_degrees[0] == 0 ? p_coeffs[0] : 0;
  else
  {
    res = 0;
    for(j = 0; j < p_num_monomials; j++)
    {
      if(p_coeffs[j])
      {
        curr_power = multiply(p_f, curr_power, power(p_f, p_y, p_degrees[j] - curr_exp));
        curr_exp = p_degrees[j];
        res ^= multiply(p_f, curr_power, p_coeffs[j]);
      }
    }
  }
  return res;
}

uint32_t evaluate_sparse_polynomial_log(
    TfiniteField* p_f,
    uint32_t* p_log_coeffs,
    uint32_t* p_degrees,
    uint32_t p_num_coeffs,
    uint32_t p_y)
{
  uint32_t j, res = 0, order = p_f->n-1;
  uint32_t logzero = p_f->logzero;
  uint32_t* exp = p_f->exp;
  if(p_y == 0)
  {
    if(p_num_coeffs > 0 && p_degrees[0] == 0 && p_log_coeffs[0] != logzero)
    {
      res = exp[p_log_coeffs[0]];
    }
  }
  else
  {
    uint64_t v;
    for(j = 0; j < p_num_coeffs; j++)
    {
      if(p_log_coeffs[j] != logzero)
      {
        // we want to compute (p_fieldPtr->log[i] * p_polyDegrees[j]) % Nfld
        // but a direct computation with uint32_t may overfow
        v = (uint64_t) p_f->log[p_y];
        v*= p_degrees[j];
        uint32_t w = (p_log_coeffs[j]+ v) % order;
        res^= exp[w];
      }
    }
  }
  return res;
}

void multiply_polynomial_by_monomial(
    TfiniteField* p_f,
    uint32_t* pio_poly,
    uint32_t p_degree,
    uint32_t p_nonzero_root_log)
{
  pio_poly[p_degree + 1] = pio_poly[p_degree];
  for (uint32_t i = p_degree; i != 0; i--)
  {
    pio_poly[i] = pio_poly[i - 1] ^ multiply_by_log(p_f, pio_poly[i], p_nonzero_root_log);
  }
  pio_poly[0] = multiply_by_log(p_f, pio_poly[0], p_nonzero_root_log);
}

void euclidean_division_sparse_reductor(
    TfiniteField* p_f,
    const uint32_t* p_reductee,
    uint32_t p_reductee_degree,
    uint32_t* p_reductor_coeffs,
    uint32_t* p_reductor_degrees,
    uint32_t p_reductor_num_terms,
    uint32_t* po_remainder,
    uint32_t* po_dividend,
    uint32_t p_dividend_factor_log)
{
  const uint32_t order = p_f->n - 1;
  const uint32_t* log = p_f->log;
  const uint32_t* exp = p_f->exp;
  uint32_t i, k;
  uint32_t reductorDegree = p_reductor_degrees[p_reductor_num_terms - 1];
  // value below may be equal to order if head coeff is 1, but the code is still correct in this case
  uint32_t reductorHeadLogInv = order - log[p_reductor_coeffs[p_reductor_num_terms - 1]];
  if (reductorDegree > p_reductee_degree)
  {
    for(i = 0; i < p_reductee_degree + 1; i++) po_remainder[i] = p_reductee[i];
    for(; i < reductorDegree; i++)             po_remainder[i] = 0;
    return;
  }
  for(i = 0; i < reductorDegree; i++) po_remainder[i] = 0;
  if (po_dividend)
  {
    for(i = 0; i < p_reductee_degree - reductorDegree + 1; i++) po_dividend[i] = 0;
  }

  uint32_t i_mod_reductor_degree = p_reductee_degree % reductorDegree;
  for (i = p_reductee_degree; i > reductorDegree - 1; i--)
  {
    // here, and in the whole loop, i_mod_reductor_degree = i mod reductorDegree;
    // introduce a new monomial of the reductee
    const uint32_t newCoeff = po_remainder[i_mod_reductor_degree] ^ p_reductee[i];
    po_remainder[i_mod_reductor_degree] = 0;
    if (newCoeff)
    {
      // add [reductor times X^(i-d)] * coeff to the reductee
      uint32_t coeffLog = log[newCoeff] + reductorHeadLogInv;
      if(coeffLog >= order) coeffLog -= order; // coeff is nonzero
      if (po_dividend)
      {
        uint32_t dividendCoeffLog = coeffLog + p_dividend_factor_log;
        if(dividendCoeffLog >= order) dividendCoeffLog -= order;
        po_dividend[i - reductorDegree] = exp[dividendCoeffLog];
      }

      for (k = 0; k < p_reductor_num_terms - 1; k++)
      {
        const uint32_t reductorCoeff = p_reductor_coeffs[k];
        // we do: po_remainderPtr[idx] += reductorCoeff * coeff
        // with idx = (p_reductorDegreesPtr[k] + i) % reductorDegree
        if(reductorCoeff)
        {
          uint32_t idx = p_reductor_degrees[k] + i_mod_reductor_degree;
          if(idx >= reductorDegree) idx -= reductorDegree;
          uint32_t sumLog = log[reductorCoeff] + coeffLog;
          if(sumLog >= order) sumLog -= order;
          po_remainder[idx] ^= exp[sumLog];
        }
      }
    }
    if(i_mod_reductor_degree == 0) i_mod_reductor_degree = reductorDegree;
    i_mod_reductor_degree--;
  }
  for (i = 0; i < reductorDegree; i++) po_remainder[i] ^= p_reductee[i];
}

void euclidean_division_sparse_reductor_ref(
    TfiniteField* p_f,
    const uint32_t* p_reductee,
    uint32_t p_reductee_degree,
    uint32_t* p_reductor_coeffs,
    uint32_t* p_reductor_degrees,
    uint32_t p_reductor_num_terms,
    uint32_t* po_remainder,
    uint32_t* po_dividend,
    uint32_t p_dividend_factor)
{
  uint32_t i, k;
  uint32_t reductorDegree = p_reductor_degrees[p_reductor_num_terms - 1];
  uint32_t reductorHead = p_reductor_coeffs[p_reductor_num_terms - 1];
  uint32_t inv = inverse(p_f, reductorHead);

  if (reductorDegree > p_reductee_degree)
  {
    for(i = 0; i < p_reductee_degree + 1; i++) po_remainder[i] = p_reductee[i];
    for(; i < reductorDegree; i++)             po_remainder[i] = 0;
    return;
  }
  for(i = 0; i < reductorDegree; i++) po_remainder[i] = 0;
  if (po_dividend)
  {
    for(i = 0; i < p_reductee_degree - reductorDegree + 1; i++) po_dividend[i] = 0;
  }

  uint32_t i_mod_reductor_degree = p_reductee_degree % reductorDegree;
  for (i = p_reductee_degree; i > reductorDegree - 1; i--)
  {
    // here, and in the whole loop, i_mod_reductor_degree = i mod reductorDegree;
    // introduce a new monomial of the reductee
    const uint32_t newCoeff = po_remainder[i_mod_reductor_degree] ^ p_reductee[i];
    po_remainder[i_mod_reductor_degree] = 0;
    if (newCoeff)
    {
      // add [reductor times X^(i-d)] * coeff to the reductee
      uint32_t coeff = multiply(p_f, newCoeff, inv);
      if (po_dividend)
      {
        po_dividend[i - reductorDegree] = multiply(p_f, coeff, p_dividend_factor);
      }

      for (k = 0; k < p_reductor_num_terms - 1; k++)
      {
        const uint32_t reductorCoeff = p_reductor_coeffs[k];
        // we do: po_remainderPtr[idx] += reductorCoeff * coeff
        // with idx = (p_reductorDegreesPtr[k] + i) % reductorDegree
        if(reductorCoeff)
        {
          uint32_t idx = p_reductor_degrees[k] + i_mod_reductor_degree;
          if(idx >= reductorDegree) idx -= reductorDegree;
          po_remainder[idx] ^= multiply(p_f, reductorCoeff, coeff);
        }
      }
    }
    if(i_mod_reductor_degree == 0) i_mod_reductor_degree = reductorDegree;
    i_mod_reductor_degree--;
  }
  for (i = 0; i < reductorDegree; i++) po_remainder[i] ^= p_reductee[i];
}
