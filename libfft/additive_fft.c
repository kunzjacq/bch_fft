#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "additive_fft.h"
#include "utils.h"

// #define MEASURE_TIME
// #define _CHECK_COMPUTATIONS

void evaluate_polynomial_additive_FFT(TfiniteField* p_f, uint32_t* p_poly, uint32_t* po_result, uint32_t p_num_terms)
{
#if 0
  additive_DFT_opt(p_f, p_poly, p_num_terms);
#else
  additive_DFT_ref(p_f, p_poly, p_num_terms);
#endif
  for (uint32_t i = 0; i < p_f->n - 1; i++) po_result[i] = p_poly[p_f->exp[i]];
}

void additive_DFT_opt(TfiniteField* p_f, uint32_t* p_poly, uint32_t p_poly_degree)
{
  uint32_t k, i, j;
  uint32_t y_log;
  const uint32_t m       = p_f->m;
  const uint32_t logzero = p_f->logzero;
  const uint32_t order   = p_f->n - 1;
  const uint32_t* log    = p_f->log;
  const uint32_t* exp    = p_f->exp;
  uint32_t* divLog  = (uint32_t*) calloc(m * (m + 1), sizeof(uint32_t));
  uint32_t* y       = (uint32_t*) calloc(m, sizeof(uint32_t));
  uint32_t* degrees = (uint32_t*) calloc(m + 1, sizeof(uint32_t));
  uint32_t* buf     = (uint32_t*) calloc(p_f->n, sizeof(uint32_t));
  for(i = 0; i < m * (m + 1); i++) divLog[i] = logzero;

#ifdef MEASURE_TIME
  double t = absolute_time();
#endif
  degrees[0] = 0;
  for(j = 1; j < m + 1; j++)
  {
    degrees[j] = 1 << (j - 1);
  }
  //divLog[i*(m+1)]         = log(constant coefficient of P_i)
  //divLog[i*(m+1) + j + 1] = log(degree 2^j coefficient of P_i)
  divLog[1] = 0; // P_0 = x^0.X
  y[0] = 1;      // x_0 = 1 => P_0(x_0) = 1

  for(i = 1; i < m; i++)
  {
    y_log = log[y[i - 1]]; // non-zero
    const uint32_t* p_i_minus_1 = divLog + (i - 1) * (m + 1);
    uint32_t* p_i = divLog + i * (m + 1);
    for(j = 1; j < i + 1; j++)
    {
      uint32_t coeff = p_i_minus_1[j];
      if(coeff == logzero) p_i[j] = logzero;
      else
      {
        coeff += y_log;
        if(coeff >= order) coeff -= order;
        p_i[j] = coeff;
      }
    }
    for(j = 2; j < m + 1; j++)
    {
      uint32_t coeff = p_i_minus_1[j - 1];
      if(coeff != logzero)
      {
        uint32_t Sqlog = coeff * 2;
        if(Sqlog >= order) Sqlog -= order;
        p_i[j] = (p_i[j] == logzero) ? Sqlog : log[exp[p_i[j]] ^ exp[Sqlog]];
      }
    }
    y[i] = evaluate_sparse_polynomial_log(p_f, p_i, degrees, i + 2, 1 << i);
  }
#ifdef MEASURE_TIME
    t = absolute_time() - t;
    printf("Polynomial computation time: %f\n", t);
    fflush(stdout);
#endif

  for(i = 0; i < m; i++)
  {
#ifdef MEASURE_TIME
    t = absolute_time();
#endif
    const uint32_t o  = 1 << (m - i);
    const uint32_t ho = 1 << (m - i - 1);
    uint32_t* remainder = buf;
    uint32_t* dividend  = buf + ho;
    // reductor = P_(m-1-i)
    uint32_t* reductorLogCoeffs = divLog + (m + 1) * (m - 1 - i);
    const uint32_t v = log[evaluate_sparse_polynomial_log(
          p_f, reductorLogCoeffs, degrees, m - i + 1, ho)];
    const uint32_t reductorDegree = degrees[m - i];
    const uint32_t reducteeDegree = p_poly_degree < o - 1 ? p_poly_degree : o - 1;
    const uint32_t reducteeDegreeReduced = reducteeDegree % reductorDegree;

    // shortcut for small degree input polynomials
    if(p_poly_degree > 0 && p_poly_degree < reductorDegree)
    {
      for(j = 0; j < (1u << i); j++)
      {
        for(k = 0; k < ho; k++) p_poly[k + j * o + ho] = p_poly[k];
      }
    }
    else
    {
      for(j = 0; j < (1u << i); j++)
      {
        uint32_t* poly_ij = p_poly + j * o;
        // change constant term of P_(m-1-i) so that P_(m-1-i)(j * o) = 0
        // first set it to zero so that the next line is evaluated on the right polynomial
        reductorLogCoeffs[0] = logzero;
        reductorLogCoeffs[0] =
            log[evaluate_sparse_polynomial_log(
              p_f, reductorLogCoeffs, degrees, m - i + 1, j * o)];

        // in the rest of iteration of the loop in j, we
        // divide poly_ij by sparse polynomial reductor.
        // if remainer is r and dividend is q,
        // poly_ij      is set to r
        // poly_ij + ho is set to r + v * q.

        // this zeroes dividend and remainder which are the two halves of buf[0 ... o-1]
        for(k = 0; k < o; k++) buf[k] = 0;
        uint32_t ii_mod_reductorDegree = reducteeDegreeReduced;
        for (uint32_t ii = reducteeDegree; ii > reductorDegree - 1; ii--)
        {
          //at all times in loop, ii_mod_reductorDegree = ii mod reductorDegree
          const uint32_t newCoeff = remainder[ii_mod_reductorDegree] ^ poly_ij[ii];
          remainder[ii_mod_reductorDegree] = 0;
          if (newCoeff)
          {
            // ratio = constant by which reductor should be multiplied to cancel reductee term =
            // reductee coeff since reductor is monic
            const uint32_t ratioLog = log[newCoeff];
            uint32_t dividendCoeffLog = ratioLog + v;
            if(dividendCoeffLog >= order) dividendCoeffLog -= order;
            dividend[ii - reductorDegree] = exp[dividendCoeffLog];
            for (k = 0; k < m - i; k++) // indexes all monomials of reductor except head
            {
              uint32_t reductorCoeffLog = reductorLogCoeffs[k];
              if(reductorCoeffLog != logzero)
              {
                uint32_t sumLog = reductorCoeffLog + ratioLog;
                if(sumLog >= order) sumLog -= order;
                uint32_t idx = ii_mod_reductorDegree + degrees[k];
                if(idx >= reductorDegree) idx -= reductorDegree;
                remainder[idx] ^= exp[sumLog];
              }
            }
          }
          //else dividend[ii - reductorDegree] = 0;
          if(ii_mod_reductorDegree == 0) ii_mod_reductorDegree = reductorDegree;
          ii_mod_reductorDegree--;
        }
        for(k = 0; k < ho; k++)
        {
          poly_ij[k]     ^= remainder[k];
          poly_ij[k + ho] = poly_ij[k] ^ dividend[k];
        }
      }
    }
#ifdef MEASURE_TIME
    t = absolute_time() - t;
    printf("Step %d time: %lf sec\n", i, t);
    fflush(stdout);
#endif
  }
  free(buf);
  free(degrees);
  free(y);
  free(divLog);
}


int check_fft_polynomials(TfiniteField* p_f, uint32_t* div, uint32_t* degrees, uint32_t* y)
{
  uint32_t i, k;
  const uint32_t m     = p_f->m;
  uint32_t* values1 = (uint32_t*) calloc(p_f->n, sizeof(uint32_t));
  uint32_t* values2 = (uint32_t*) calloc(p_f->n, sizeof(uint32_t));
  printf("P_0: \n");
  for(k = 0; k < p_f->n; k++)
  {
    values1[k] = evaluate_sparse_polynomial(p_f, div, degrees, 2, k);
  }
  int ok = 1, global_ok = 1;
  if(values1[0] != 0) ok = 0;
  for(k = 1; k < (1u << m); k++) if(values1[k] == 0) ok = 0;
  if(!ok) printf("Polynomial P_0 has wrong roots\n");
  for(i = 1; i < m; i++)
  {
#if 0
    printf("y_%d=0x%02x\n", i, y[i]);
    printf("P_%02d: ", i);
    int has_coeff = 0;
    for(k = 0; k < i + 2; k++)
    {
      if(has_coeff) printf(" + ");
      uint32_t coeff = div[(m + 1) * i + k];
      if(coeff)
      {
        has_coeff = 1;
        if(coeff != 1) printf("0x%02x.X^%d ", div[(m + 1) * i + k], degrees[k]);
        else printf("X^%d ", degrees[k]);
      }
    }
#endif
    printf("\n");
    ok = 1;
    int disc = 0;

    uint32_t* t = values2;
    values2 = values1;
    values1 = t;
    for(k = 0; k < p_f->n; k++)
    {
      values1[k] = evaluate_sparse_polynomial(p_f,  div + (m + 1) * i, degrees, i + 2, k);
    }

    for(k = 0; k < (1u << m); k++)
    {
      uint32_t d = multiply(p_f, values2[k], values2[k^(1<<(i-1))]) != values1[k];
      if(d)
      {
        ok = 0;
        disc++;
      }
    }
    if(ok) printf("ok: P_%d = P_%d(X) x P_%d(X-x_%d)\n", i, i - 1, i - 1, i - 1);
    else   printf("nok: P_%d != P_%d(X) x P_%d(X-x_%d) (%d discrepancies)\n", i, i - 1, i - 1, i - 1, disc);
    global_ok &= ok;
    ok = 1;
    for(k = 0; k < (1u << i); k++)
    {
      if(values1[k] != 0)
      {
        printf("%d+ ",k);
        ok = 0;
      }
    }
    for(; k < (1u << m); k++)
    {
      if(values1[k] == 0)
      {
        printf("%d- ",k); ok = 0;
      }
    }
    printf("\n");
    if(ok) printf("ok: Polynomial P_%d has the expected roots\n",i);
    else   printf("nok: Polynomial P_%d has wrong roots\n",i);
    global_ok &= ok;
  }
  free(values2);
  free(values1);
  return global_ok;
}

void additive_DFT(TfiniteField* p_f, uint32_t* p_poly, uint32_t p_poly_degree)
{
  const uint32_t m     = p_f->m;
  const uint32_t order = p_f->n - 1;
  const uint32_t* log  = p_f->log;
  const uint32_t* exp  = p_f->exp;
  uint32_t* div      = (uint32_t*) calloc(m * (m + 1), sizeof(uint32_t));
  uint32_t* y        = (uint32_t*) calloc(m,           sizeof(uint32_t));
  uint32_t* degrees  = (uint32_t*) calloc(m + 1,       sizeof(uint32_t));
  uint32_t* buf      = (uint32_t*) calloc(p_f->n,      sizeof(uint32_t));
  uint32_t i, j, k, ylog;

#ifdef MEASURE_TIME
    double t = absolute_time();
#endif
  degrees[0] = 0;
  for(j = 1; j < m + 1; j++) degrees[j] = 1 << (j - 1);
  // div[i * (m + 1)]         = constant coefficient of P_i
  // div[i * (m + 1) + j + 1] = degree 2**j coefficient of P_i
  div[1] = 1; // P_0 = X
  y[0]   = 1; // x_0 = 1 => P_0(x_0) = 1

  // recursive computation of the P_i which are linear polynomials
  // P_i vanishes on 0 ... 2**i - 1 (included)
  // y_ i = P_i(x_i), x_i = 2**i
  // P_i+1(X) = P_i(X) * P_i(X - x_i)
  //          = P_i(X) * (P_i(X) - y_i)
  //          = P_i**2(X) - y_i P_i(X)
  for(i = 1; i < m; i++)
  {
    const uint32_t* p_i_minus_1 = div + (i - 1) * (m + 1);
    uint32_t*               p_i = div + i       * (m + 1);
    ylog = log[y[i - 1]]; // non-zero
    // add y_i P_i-1
    for(j = 1; j < i + 1; j++)
    {
      const uint32_t coeff = p_i_minus_1[j];
      p_i[j] = multiply_by_log(p_f, coeff, ylog);
    }
    // add P_i**2
    for(j = 2; j < m + 1; j++)
    {
      // 2 * degree of  j-1 th monomial of P_i-1 = degree of j-th monomial of P_i
      const uint32_t coeff = p_i_minus_1[j - 1];
      if(coeff != 0)
      {
        uint32_t Sqlog = log[coeff] * 2;
        if(Sqlog >= order) Sqlog -= order;
        p_i[j] ^= exp[Sqlog];
      }
    }
    // y_i = p_i (1 << i)
    y[i] = evaluate_sparse_polynomial(p_f, p_i, degrees, i + 2, 1 << i);
  }
#ifdef MEASURE_TIME
    t = absolute_time() - t;
    printf("Polynomial computation time (could be precomputed): %lf\n", t);
#endif

#ifdef _CHECK_COMPUTATIONS
    check_fft_polynomials(p_f, div, degrees, y);
#endif

  for(i = 0; i < m; i++)
  {
#ifdef MEASURE_TIME
    t = absolute_time();
#endif
    const uint32_t o  = 1 << (m - i);
    const uint32_t ho = 1 << (m - i - 1);
    //reductor = P_(m - 1 - i)
    uint32_t* reductorCoeffs = div + (m + 1) * (m - 1 - i);
    // v = reductor(ho)
    const uint32_t reducteeDegree = p_poly_degree < o - 1 ? p_poly_degree : o - 1;
    const uint32_t reductorDegree = degrees[m - i];
    const uint32_t log_v =
        log[evaluate_sparse_polynomial(p_f, reductorCoeffs, degrees, m - i + 1, ho)];
    // shortcut for small degree input polynomials:
    // if polynomial degree is smaller than reductor degree, no reduction occurs
    if(p_poly_degree > 0 && p_poly_degree < reductorDegree)
    {
      for(j = 0; j < (1u << i); j++)
      {
        for(k = 0; k < ho; k++) p_poly[k + j * o + ho] = p_poly[k];
      }
    }
    else
    {
      for(j = 0; j < (1u << i); j++)
      {
        uint32_t* poly_ij = p_poly + j * o;
        //change constant term of P_(m-1-i) so that P_(m-1-i)(j * o) = 0
        reductorCoeffs[0] ^=
            evaluate_sparse_polynomial(
              p_f, reductorCoeffs, degrees, m - i + 1, j * o);

        // divide poly_ij by sparse polynomial reductor.
        // buf will hold remainer r
        // buf + ho will hold (q x v) where q is the dividend.
        euclidean_division_sparse_reductor(
              p_f, poly_ij, reducteeDegree, reductorCoeffs, degrees,
              m - i + 1, buf, buf + ho, log_v);

        for(k = 0; k < ho; k++)
        {
          poly_ij[k]      = buf[k];                   // r
          poly_ij[k + ho] = poly_ij[k] ^ buf[k + ho]; // r + v * q
        }
      }
    }
#ifdef MEASURE_TIME
    t = absolute_time() - t;
    printf("Step %d time: %lf sec\n", i, t);
#endif
  }
  free(buf);
  free(degrees);
  free(y);
  free(div);
}


void additive_DFT_ref(TfiniteField* p_f, uint32_t* p_poly, uint32_t p_poly_degree)
{
  const uint32_t m     = p_f->m;
  uint32_t* div      = (uint32_t*) calloc(m * (m + 1), sizeof(uint32_t));
  uint32_t* y        = (uint32_t*) calloc(m,           sizeof(uint32_t));
  uint32_t* degrees  = (uint32_t*) calloc(m + 1,       sizeof(uint32_t));
  uint32_t* buf      = (uint32_t*) calloc(p_f->n,      sizeof(uint32_t));
  uint32_t i, j, k, yval;

  degrees[0] = 0;
  for(j = 1; j < m + 1; j++) degrees[j] = 1 << (j - 1);
  // div[i * (m + 1)]         = constant coefficient of P_i
  // div[i * (m + 1) + j + 1] = degree 2**j coefficient of P_i
  div[1] = 1; // P_0 = X
  y[0]   = 1; // x_0 = 1 => P_0(x_0) = 1

  // recursive computation of the P_i which are linear polynomials
  // P_i vanishes on 0 ... 2**i - 1 (included)
  // y_ i = P_i(x_i), x_i = 2**i
  // P_i+1(X) = P_i(X) * P_i(X - x_i)
  //          = P_i(X) * (P_i(X) - y_i)
  //          = P_i**2(X) - y_i P_i(X)
  for(i = 1; i < m; i++)
  {
    const uint32_t* p_i_minus_1 = div + (i - 1) * (m + 1);
    uint32_t*               p_i = div + i       * (m + 1);
    yval = y[i - 1]; // non-zero
    // add y_i P_i-1
    for(j = 1; j < i + 1; j++)
    {
      const uint32_t coeff = p_i_minus_1[j];
      p_i[j] = multiply(p_f, coeff, yval);
    }
    // add P_i**2
    for(j = 2; j < m + 1; j++)
    {
      // 2 * degree of  j-1 th monomial of P_i-1 = degree of j-th monomial of P_i
      const uint32_t coeff = p_i_minus_1[j - 1];
      p_i[j] ^= square(p_f, coeff);
    }
    // y_i = p_i (1 << i)
    y[i] = evaluate_sparse_polynomial_ref(p_f, p_i, degrees, i + 2, 1 << i);
  }

#ifdef _CHECK_COMPUTATIONS
  check_fft_polynomials(p_f, div, degrees, y);
#endif
  for(i = 0; i < m; i++)
  {
    const uint32_t o  = 1 << (m - i);
    const uint32_t ho = 1 << (m - i - 1);
    //reductor = P_(m - 1 - i)
    uint32_t* reductorCoeffs = div + (m + 1) * (m - 1 - i);
    // v = reductor(ho)
    const uint32_t reducteeDegree = p_poly_degree < o - 1 ? p_poly_degree : o - 1;
    const uint32_t reductorDegree = degrees[m - i];
    const uint32_t v =
        evaluate_sparse_polynomial_ref(p_f, reductorCoeffs, degrees, m - i + 1, ho);
    // shortcut for small degree input polynomials:
    // if polynomial degree is smaller than reductor degree, no reduction occurs
    if(p_poly_degree > 0 && p_poly_degree < reductorDegree)
    {
      for(j = 0; j < (1u << i); j++)
      {
        for(k = 0; k < ho; k++) p_poly[k + j * o + ho] = p_poly[k];
      }
    }
    else
    {
      for(j = 0; j < (1u << i); j++)
      {
        uint32_t* poly_ij = p_poly + j * o;
        //change constant term of P_(m-1-i) so that P_(m-1-i)(j * o) = 0
        reductorCoeffs[0] ^= evaluate_sparse_polynomial_ref(p_f, reductorCoeffs, degrees, m - i + 1, j * o);

        // divide poly_ij by sparse polynomial reductor.
        // buf will hold remainer r
        // buf + ho will hold (q x v) where q is the dividend.
        //euclidean_division_sparse_reductor(
        //      p_f, poly_ij, reducteeDegree, reductorCoeffs, degrees,
        //      m - i + 1, buf, buf + ho, p_f->log[v]);
        euclidean_division_sparse_reductor_ref(
              p_f, poly_ij, reducteeDegree, reductorCoeffs, degrees,
              m - i + 1, buf, buf + ho, v);

        for(k = 0; k < ho; k++)
        {
          poly_ij[k]      = buf[k];                   // r
          poly_ij[k + ho] = poly_ij[k] ^ buf[k + ho]; // r + v * q
        }
      }
    }
  }

  free(buf);
  free(degrees);
  free(y);
  free(div);
}

