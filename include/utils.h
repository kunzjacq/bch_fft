#include <stdint.h>
#include <stddef.h>

#if defined(__cplusplus)
#define UNUSED(x)
#elif defined(__GNUC__)
#define UNUSED(x)       x __attribute__((unused))
#else
#define UNUSED(x)       x
#endif

double absolute_time();
void init_time();
uint32_t urand();
int compare_results(uint32_t* p_dataA, uint32_t* p_dataB, uint32_t p_groupSize, uint32_t szLog);

