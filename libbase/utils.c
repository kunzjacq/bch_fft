#include "utils.h"

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>

static long start_time_sec = 0;

void init_time()
{
  struct timeval start;
  gettimeofday(&start, NULL);
  start_time_sec = start.tv_sec;
}

double absolute_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  // subtracting start_time_sec aims at improving the precision of the conversion to double
  return (t.tv_sec - start_time_sec) * 1. + t.tv_usec *1e-6;
}

uint32_t urand(){
  return (uint32_t) rand();
}

int compare_results(uint32_t* p_dataA, uint32_t* p_dataB, uint32_t p_groupSize, uint32_t szLog)
{
  uint32_t i, checksumA, checksumB;
  int error = 0;
  checksumA = 0;
  checksumB = 0;
  for (i = 0; i < p_groupSize; i++)
  {
    checksumB += p_dataB[i];
    checksumA += p_dataA[i];
    if (p_dataA[i] != p_dataB[i]) error++;
  }

  char szStr[] = "%04x ";
  if(4 <= szLog && szLog <= 32) {
    szStr[2]='0'+((char)szLog + 3) / 4;
  }

  if(error)
  {
    printf("************* %d errors ********************\n", error);
    if(p_groupSize > 20)
    {
      uint32_t numdisplay = 10;
      for(i = 0; i < numdisplay; i++) printf(szStr, p_dataA[i]);
      printf("... ");
      for(i = 0; i < numdisplay; i++) printf(szStr, p_dataA[p_groupSize - numdisplay + i]);
      printf("\n");
      for(i = 0; i < numdisplay; i++) printf(szStr, p_dataB[i]);
      printf("... ");
      for(i = 0; i < numdisplay; i++) printf(szStr, p_dataB[p_groupSize - numdisplay + i]);
      printf("\n");
    }
    else
    {
      for(i = 0; i < p_groupSize; i++) printf("%04x ", p_dataA[i]);
      printf("\n");
      for(i = 0; i < p_groupSize; i++) printf("%04x ", p_dataB[i]);
      printf("\n");
    }
    printf("checksums: %8x / %8x\n", checksumA, checksumB);
  }
  return error;
}
