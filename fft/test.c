#include "finite_field.h"
#include "additive_fft.h"
#include "multiplicative_fft.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* log2 of finite field sizes for which we will not bench multiplicative
 * fft timing since the multiplicative group size is not smooth and therefore
 * fft performance is poor */
static const uint32_t blacklist[] = {13, 17, 19, 23, 26};
static const uint32_t blacklist_size = 5;
static const uint32_t minSize = 14;
static const uint32_t maxSize = 21;
// plain fourier transform is slow, so don't explore large sizes
static const uint32_t maxSize_plain_FT = 14;

int partial_dft_one_test(TfiniteField *f,  uint32_t* p_buffers, uint32_t numCoeffs);
int multiplicative_dft_test(uint32_t* p_buffers);
int additive_dft_test(uint32_t* p_buffers);
int partial_dft_test(uint32_t* p_buffers);

int main(int UNUSED(argc), char* UNUSED(argv[]))
{
  uint32_t* buffers = (uint32_t*) calloc(6 * (1 << maxSize), 4);
  int error = 0;
  init_time();

  error |= multiplicative_dft_test(buffers);
  error |= additive_dft_test(buffers);
  error |= partial_dft_test(buffers);

  free(buffers);
  if (!error)
  {
    printf("\nTest successful!\n");
    return EXIT_SUCCESS;
  }
  else
  {
    printf("\n********** Test failed! **********\n");
    return EXIT_FAILURE;
  }
}

int multiplicative_dft_test(uint32_t* p_buffers)
{
  uint32_t i,szLog, blacklist_idx = 0;
  int error = 0, local_error;
  double t1, t2;
  uint32_t* refIn, *refOut, *buffer1, *buffer2;
  refIn = p_buffers;
  for(szLog = minSize; szLog <= maxSize_plain_FT; szLog++)
  {
    while(blacklist_idx < blacklist_size && szLog > blacklist[blacklist_idx]) blacklist_idx++;
    if(blacklist_idx < blacklist_size && szLog == blacklist[blacklist_idx]) continue;
    TfiniteField f = create_binary_finite_field(szLog);
    uint32_t repeat = f.n <= (1u << maxSize_plain_FT) ?  (1u << maxSize_plain_FT) / f.n : 1;
    uint32_t mult_order = f.n - 1;
    refOut  = p_buffers +     f.n;
    buffer1 = p_buffers + 2 * f.n;
    buffer2 = p_buffers + 3 * f.n;
    printf("Computing a FFT of size 2^%d - 1 = %d\n", szLog, mult_order);
    fflush(stdout);
    srand(5);
    for (i = 0; i < mult_order; i++) refIn[i] = i;

    printf("Performing Fast Fourier Transform using the factoring of %d\n", mult_order);
    t2 = absolute_time();
    for(i = 0 ; i < repeat; i++)
    {
      memcpy(refOut, refIn, f.n * 4);
      multiplicative_FFT(&f, refOut, buffer1);
    }
    t2 = absolute_time() - t2;
    printf("time: %.2f sec.\n", t2);

    printf("Computing Fourier transform with no speedup\n");
    t1 = absolute_time();
    for(i = 0 ; i < repeat; i++)
    {
      memcpy(buffer2, refIn, f.n * 4);
      DFT_direct(&f, buffer2, buffer1);
    }
    t1 = absolute_time() - t1;
    local_error = compare_results(buffer2, refOut, mult_order, szLog);
    printf("time: %.2f sec.\n", t1);
    printf("speedup : %.2f\n", t1/t2);
    if(local_error) printf("********** Multiplicative DFT: wrong result **********\n");
    else printf("Success!\n");
    error |= local_error;
    free_finite_field(&f);
  }
  if(error) printf("********** Multiplicative DFT Failed **********\n");
  return error;
}

int additive_dft_test(uint32_t* p_buffers)
{
  uint32_t i;
  uint32_t szLog[] = {8, 16, 20};
  double t1, t2;
  int local_error, error = 0;
  srand(5);
  for(unsigned int j = 0; j < sizeof(szLog)/4; j++)
  {
    if(szLog[j] < minSize || szLog[j] > maxSize) continue;
    TfiniteField f = create_binary_finite_field(szLog[j]);
    uint32_t mult_order = f.n - 1;
    printf("\nTest of additive DFT of size %d\n", mult_order);
    uint32_t* refIn, *refOut, *buffer1, *buffer2;
    refIn = p_buffers;
    refOut  = p_buffers +     f.n;
    buffer1 = p_buffers + 2 * f.n;
    buffer2 = p_buffers + 3 * f.n;

    uint32_t max_degree = f.n - 1;

    for (i = 0; i < f.n; i++)
    {
      refIn[i] = i <= max_degree ? i & 0x3FF : 0;
      //refIn[i] = urand() & mult_order;
    }
    refIn[f.n - 1] = 0;
    uint32_t repeat = 1;//f.n <= (1u<<20) ?  (1u<<20) / f.n : 1;
    printf("Doing full DFT...\n");
    t1 = absolute_time();
    for(i = 0 ; i < repeat; i++) {
      memcpy(refOut, refIn, f.n * 4);
      multiplicative_FFT(&f, refOut, buffer1);
    }
    t1 = absolute_time() - t1;
    printf("time of full DFT: %.2f sec.\n", t1);
    printf("Doing additive DFT...\n");
    t2 = absolute_time();
    for(i = 0 ; i < repeat; i++)
    {
      memcpy(buffer1, refIn, f.n * 4);
      evaluate_polynomial_additive_FFT(&f, buffer1, buffer2, max_degree);
    }
    t2 = absolute_time() - t2;
    printf("time of additive DFT: %.2f sec.\n", t2);
    printf("speed ratio : %.2f\n", t1/t2);
    fflush(stdout);
    local_error = compare_results(buffer2, refOut, mult_order, szLog[j]);
    error |= local_error;
    if(local_error) printf("********** Additive DFT: wrong result **********\n");
    else printf("Success!\n");
    free_finite_field(&f);
  }
  if(error)
  {
    printf("********** Additive DFT Failed **********\n");
  }
  else
  {
    printf("additive DFT succeeded\n");
  }
  return error;
}

int partial_dft_test(uint32_t* p_buffers)
{
  uint32_t szLog, blacklist_idx = 0, mult_order;
  int local_error = 0;
  printf("\nTest of small degree polynomial evaluation with partial DFT\n");
  for(szLog = minSize; szLog <= maxSize; szLog++)
  {
    while(blacklist_idx < blacklist_size && szLog > blacklist[blacklist_idx]) blacklist_idx++;
    if(blacklist_idx < blacklist_size && szLog == blacklist[blacklist_idx]) continue;
    printf("Field size: 2^%d\n", szLog);
    TfiniteField f = create_binary_finite_field(szLog);
    mult_order = f.n - 1;
    uint32_t lastfactor = factors(szLog)[numFactors(szLog)-1];
    uint32_t numCoeffs = mult_order/lastfactor;
    srand(5);
    local_error |= partial_dft_one_test(&f, p_buffers, numCoeffs);
    if(numFactors(szLog) >= 2)
    {
      uint32_t prevfactor = factors(szLog)[numFactors(szLog) - 2];
      numCoeffs = mult_order/(lastfactor*prevfactor);
      local_error |= partial_dft_one_test(&f, p_buffers, numCoeffs);
    }

    free_finite_field(&f);
    printf("\n");
  }
  if(local_error)
  {
    printf("********** small degree polynomial evaluation with partial DFT Failed **********\n");
  }
  else
  {
    printf("small degree polynomial evaluation with partial DFT succeeded\n");
  }
  return local_error;
}

int partial_dft_one_test(TfiniteField *f,  uint32_t* p_buffers, uint32_t numCoeffs)
{
  uint32_t mult_order = f->n - 1;
  uint32_t i;
  double t1, t2;
  int local_error;
  uint32_t* refIn, *refOut, *buffer1, *buffer2, * result;
  refIn   = p_buffers;
  refOut  = p_buffers +     f->n;
  buffer1 = p_buffers + 2 * f->n;
  buffer2 = p_buffers + 3 * f->n;
  result  = p_buffers + 4 * f->n;

  for (i = 0; i < f->n; i++)
  {
    refIn[i]  = 0;
    refOut[i] = 0;
  }
  for (i = 0; i < numCoeffs; i++)
  {
    refIn[i] = ((unsigned int) rand()) & mult_order;
  }
  printf("Polynomial degree: %d\n", numCoeffs - 1);
  uint32_t repeat = f->n <= (1u<<20) ?  (1u<<20) / f->n : 1;
  printf("Doing full DFT...\n");
  t1 = absolute_time();
  for(i = 0 ; i < repeat; i++) {
    memcpy(buffer1, refIn, f->n * 4);
    multiplicative_FFT(f, buffer1, buffer2);
  }
  memcpy(refOut, buffer1, mult_order * 4);
  t1 = absolute_time() - t1;
  printf("time of full DFT: %.2f sec.\n", t1);

  printf("Computing a partial DFT with evaluate_polynomial_multiplicative_FFT()...\n");
  t2 = absolute_time();
  for(i = 0 ; i < repeat; i++)
  {
    evaluate_polynomial_multiplicative_FFT(f, refIn, result, numCoeffs, buffer1);
  }
  t2 = absolute_time() - t2;
  printf("time of polynomial evaluation with evaluatePolynomial: %.2f sec.\n", t2);
  printf("speedup : %.2f\n", t1 / t2);
  local_error = compare_results(result, refOut, mult_order, f->m);
  if(!local_error) printf("Success!\n");
  fflush(stdout);
  return local_error;
}
