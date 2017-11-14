/*	Cpt S 411, Introduction to Parallel Computing
 *	School of Electrical Engineering and Computer Science
 *
 *  Joshua Shepherd
 *  2017-11-13
 * */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
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

int main(int argc, char** argv, char** envp) {

  if (argc < 3) {
    printf("Usage: %s [# of iterations] [# of threads]\n", argv[0]);
  }

  srand(0);

  uint32_t rank = 0;
  uint64_t num_rands = 1024;
  uint64_t rands_in_circ = 0;
  const uint64_t circ_radius_sq = (1 << 31) * (1 << 31);
  for (uint64_t i_iter=0; i_iter < num_rands; ++i_iter) {
    uint64_t x_rand = rand_r(&rank);
    uint64_t y_rand = rand_r(&rank);
    uint64_t dist_2_orig_sq = x_rand * x_rand + y_rand * y_rand;
    if (dist_2_orig_sq <= circ_radius_sq) {
      ++rands_in_circ;
    }
  }

  double exp_pi_ratio = rands_in_circ / num_rands;

  printf("experimental ratio: %lu / %lu\n", rands_in_circ, num_rands);
  printf("experimental ratio: %lf\n", exp_pi_ratio);

  
  return 0;
}
