/*	Cpt S 411, Introduction to Parallel Computing
 *	School of Electrical Engineering and Computer Science
 *
 *  Joshua Shepherd
 *  2017-11-13
 * */


#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include <assert.h>
#include <sys/time.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#include <algorithm>
#include <string>

uint64_t bal_rand_r(uint32_t* seed) {
  int randval = rand_r(seed);
  while (randval == (-1 << 31)) {
    randval = rand_r(seed);
  }
  return randval;
}

int main(int argc, char** argv, char** envp) {

  if (argc < 3) {
    printf("Usage: %s [# of iterations] [# of threads]\n", argv[0]);
  }

  srand(0);

  uint32_t rank = 0;
  uint64_t num_rands = std::stoull(argv[1]);//1024;
  uint64_t rands_in_circ = 0;
  const uint64_t circ_radius_sq = (((uint64_t) 1 << 31) - 1) * (((uint64_t) 1 << 31) - 1);
  uint64_t step = 100000000; // 100 million
  uint64_t i_iter = 0;
  uint32_t istep = 0;
  printf("iterations = %lu ; step = %lu\n", num_rands, step);
  while (i_iter < num_rands) {
    uint64_t step_stop = std::min(i_iter + step, num_rands);
    for (; i_iter < step_stop; ++i_iter) {
      uint64_t x_rand = bal_rand_r(&rank);
      uint64_t y_rand = bal_rand_r(&rank);
      uint64_t dist_2_orig_sq = x_rand * x_rand + y_rand * y_rand;
      if (dist_2_orig_sq <= circ_radius_sq) { // do we count the border? TODO
        ++rands_in_circ;
      }
    }
    double exp_pi_ratio = (double) 4 * rands_in_circ / step_stop;
    printf("step %u\n", istep);
    printf("  experimental ratio: 4 * %lu / %lu\n", rands_in_circ, step_stop);
    printf("  experimental ratio: %.10lf\n", exp_pi_ratio);
    ++istep;
  }

  // best case
  /*

  |xxxx
  |    xx
  |     xx
  |       x
  |_______x___

  */

  const uint64_t num_in_quadrant = (((uint64_t) 1 << 31) - 1) * (((uint64_t) 1 << 31) - 1);
  rands_in_circ = 0;
  uint64_t x = 0;
  
  // for 2^31 - 1 > y >= 0 
  //https://stackoverflow.com/questions/9044059/unsigned-integers-in-c-for-loops
  for (uint64_t y = ((uint64_t) 1 << 31) - 1; y-- != 0; ) {
    for (; x <= ((uint64_t) 1 << 31) - 1; ++x) {
      uint64_t dist_2_orig_sq = x*x + y*y;
      if(dist_2_orig_sq <= circ_radius_sq) {
        ++rands_in_circ;
        continue;
      } else {
        break;
      }
    }
    rands_in_circ += x;
  }
  
  double exp_pi_ratio = (double) 4 * rands_in_circ / num_in_quadrant;
  printf("ideal\n");
  printf("  ratio: 4 * %lu / %lu\n", rands_in_circ, num_in_quadrant);
  printf("  ratio: %.10lf\n", exp_pi_ratio);

  return 0;
}
