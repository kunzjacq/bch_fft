#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "gf2_polynomials.h"
#include "gf2_polynomials_misc.h"
#include "bch.h"
#include "utils.h"

#define BUFFER_LENGTH 512

__attribute__((__noreturn__)) void usage(char* argv[]);
int test_c(unsigned int p_length, unsigned int p_t, unsigned int p_numerr);

int main(int argc, char* argv[])
{
  unsigned int length, t, numerr;

  // ignore first argument
  // argc--; argv++;

  if (argc != 4) usage(argv);
  else
  {
    length = (unsigned int) atoi(argv[1]);
    t      = (unsigned int) atoi(argv[2]);
    numerr = (unsigned int) atoi(argv[3]);
  }
  if (t == 0 || length <= t || length > (1u << g_maxDegree)) usage(argv);
  int success = test_c(length, t, numerr);
  if (success)
  {
    printf("Test successful\n");
    return EXIT_SUCCESS;
  }
  else
  {
    printf("Test failed");
    return EXIT_FAILURE;
  }
}

int test_c(unsigned int p_length, unsigned int p_t, unsigned int p_numerr)
{
  int success = 1;
  double generationTime, encodingTime, syndromeComputationTime, elpComputationTime,
      errorLocationTime = 0., decodingTime;
  unsigned int d, tmp;
  unsigned int seed, decerror;
  uint32_t* data;
  unsigned long* encodedData, * encodedDataRef;
  uint32_t num_corrected_errors = 0, i, data_size;
  TbchCode code;
  printf("Test program with C interface\n\n");
  encodedData = 0;
  encodedDataRef = 0;
  data = 0;
  d = 2 * p_t + 1;
  // Compute the generator polynomial of BCH code
  init_time();
  generationTime = absolute_time();
  code = bch_gen_code(p_length, d);
  generationTime = absolute_time() - generationTime;
  if (code.m_packed_gen_poly == NULL)
  {
    printf("Code could not be created: codeword too small given the requested error-correcting "
        "capability");
    exit(2);
  }
  printf("Created a BCH code\n with block size %d\n with %d redundancy bits\n able to correct %d errors\n",
         code.m_length , code.m_gen_poly_degree, p_t);
  // generate random data
  seed = 1;
  srand(seed);
  data_size = code.m_length - code.m_gen_poly_degree;
  data = (uint32_t*) malloc(data_size * sizeof (int));
  for (i = 0; i < data_size; i++) data[i] = rand() & 1;
  printf("Encoding data\n");
  size_t polyBytesize = degree_to_packed_size(code.m_length - 1) * sizeof (unsigned long);
  encodedData = (unsigned long*) calloc(polyBytesize, 1);
  encodedDataRef = (unsigned long*) calloc(polyBytesize, 1);
  encodingTime = absolute_time();
  bch_encode(&code, encodedData, data, data_size);
  encodingTime = absolute_time() - encodingTime;
  memcpy(encodedDataRef, encodedData, polyBytesize);
  printf("Adding %d errors\n", p_numerr);
  i = 0;
  while(i < p_numerr)
  {
    tmp = urand() % code.m_length;
    if(get_bit(encodedDataRef, tmp) == get_bit(encodedData, tmp))
    {
      add_bit(encodedData, tmp);
      i++;
    }
  }
  printf("Decoding...\n");
  int32_t syn_error;
  uint32_t deg;
  TbchCodeBuffer buffer = bch_alloc_buffers(&code);
  printf(" Computing syndromes...\n");
  syndromeComputationTime = absolute_time();
  syn_error = bch_eval_syndrome(&code, &buffer, encodedData);
  // values of codeword on roots of generating polynomial output in buffer.m_syndrome
  syndromeComputationTime = absolute_time() - syndromeComputationTime;
  if (syn_error)
  {
    printf(" Computing error location polynomial (ELP)...\n");
    elpComputationTime = absolute_time();
    deg = bch_compute_elp(&code, &buffer);
    elpComputationTime = absolute_time() - elpComputationTime;
    printf(" ELP has degree %d\n", deg);

    if (deg <= (code.m_distance - 1) / 2)
    {
      // Find the roots of the elp that are in the index range
      // 0 ... p_codePtr->m_length - 1
      errorLocationTime = absolute_time();
      printf(" Locating and correcting errors...\n");
      num_corrected_errors = bch_locate_errors(&code, &buffer, deg);
      errorLocationTime = absolute_time() - errorLocationTime;
      for (i = 0; i < num_corrected_errors; i++) add_bit(encodedData, buffer.m_error[i]);
      if (num_corrected_errors != deg)
      {
        printf(" BCH: ELP has degree %d, but found %d roots\n", deg, num_corrected_errors);
      }
    }
    else
    {
      printf(" BCH: Too many errors, can't correct\n");
    }
  }
  else
  {
    elpComputationTime = 0;
    errorLocationTime = 0;
  }
  decodingTime = syndromeComputationTime + elpComputationTime + errorLocationTime;
  bch_free_buffers(&buffer);
  if (num_corrected_errors)
  {
    printf("According to the decoder, %d errors were corrected\n", num_corrected_errors);
  }
  else
  {
    printf("According to the decoder, some errors could not be corrected\n");
  }
  printf("Checking decoded word\n");
  decerror = 0;

  for (i = 0; i < code.m_length; i++)
  {
    if (get_bit(encodedDataRef, i) != get_bit(encodedData, i))
    {
      decerror++;
    }
  }
  if (decerror)
  {
    printf("There were %d decoding errors\n", decerror);
    success = 0;
  }
  else
  {
    printf("There are no remaining errors\n");
  }
  printf("\nGeneration time: %.2fs.\n", generationTime);
  printf("Encoding time:   %.2fs.\n", encodingTime);
  printf("Decoding time:   %.2fs.\n", decodingTime);
  printf(" syndrome evaluation time:  %.2fs.\n", syndromeComputationTime);
  printf(" ELP computation time:      %.2fs.\n", elpComputationTime);
  printf(" error location time:       %.2fs.\n", errorLocationTime);
  free(data);
  free(encodedData);
  free(encodedDataRef);
  bch_free_code(&code);
  return success;
}

__attribute__((__noreturn__)) void usage(char* argv[])
{
  printf("Usage: %s block_size error_correcting_capability number_of_errors\n", argv[0]);
  uint32_t maxValue = 1u << g_maxDegree;
  printf(" Maximum value for block size is %d\n", maxValue);
  exit(1);
}
