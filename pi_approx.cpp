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

#include <omp.h>


typedef struct rational {
  uint64_t num;
  uint64_t den;
} Rational;

uint64_t bal_rand_r(uint32_t* seed) {
  int32_t randval = rand_r(seed);
  while (randval == (-1 << 31)) {
    randval = rand_r(seed);
  }
  return randval;
}

Rational calc_ideal_ratio() {
  const uint64_t circ_radius_sq  = (((uint64_t) 1 << 31) - 1) * (((uint64_t) 1 << 31) - 1);
  const uint64_t num_in_quadrant = (((uint64_t) 1 << 31) - 1) * (((uint64_t) 1 << 31) - 1);
  uint64_t rands_in_circ = 0;
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


  Rational ideal_pi;
  ideal_pi.num = 4 * rands_in_circ;
  ideal_pi.den = num_in_quadrant;

  return ideal_pi;
}

Rational calc_exper_ratio(uint32_t rank, uint64_t num_rands) {
  uint32_t seed = rank;
  uint64_t rands_in_circ = 0;
  const uint64_t circ_radius_sq = (((uint64_t) 1 << 31) - 1) * (((uint64_t) 1 << 31) - 1);
  uint64_t step = 100000000; // 100 million
  uint64_t i_iter = 0;
  uint32_t istep = 0;
  printf("iterations = %lu ; step = %lu\n", num_rands, step);
  while (i_iter < num_rands) {
    uint64_t step_stop = std::min(i_iter + step, num_rands);
    for (; i_iter < step_stop; ++i_iter) {
      uint64_t x_rand = bal_rand_r(&seed);
      uint64_t y_rand = bal_rand_r(&seed);
      uint64_t dist_2_orig_sq = x_rand * x_rand + y_rand * y_rand;
      if (dist_2_orig_sq <= circ_radius_sq) { // do we count the border? TODO
        ++rands_in_circ;
      }
    }
    double exp_pi_ratio = (double) 4 * rands_in_circ / step_stop;

    printf("r%2u step %u: %.10lf ~= 4 * %lu / %lu \n", rank, istep, exp_pi_ratio, rands_in_circ, step_stop);

    #if 0
    printf("r%2d step %u\n", rank, istep);
    printf("r%2d  experimental ratio: 4 * %lu / %lu\n", rank, rands_in_circ, step_stop);
    printf("r%2d  experimental ratio: %.10lf\n", rank, exp_pi_ratio);
    #endif
    ++istep;
  }

  Rational approx_pi;
  approx_pi.num = 4 * rands_in_circ;
  approx_pi.den = num_rands;

  return approx_pi;
}

int main(int argc, char** argv, char** envp) {

  if (argc != 3) {
    printf("Usage: %s [# of iterations] [# of threads]\n", argv[0]);
    return -1;
  }
  
  int num_threads = atoi(argv[2]);
  #pragma omp parallel num_threads(num_threads)
  assert( num_threads == omp_get_num_threads());
  printf("# of threads: %d\n", num_threads);

  srand(0);

  //uint32_t rank = 0;
  uint64_t num_rands = std::stoull(argv[1]);//1024;


  Rational sum_approx_pi = {0,0};
  #pragma omp parallel for
  for (int tid=0; tid < num_threads; ++tid) {
    Rational indiv_approx_pi = calc_exper_ratio(tid, num_rands);
    #pragma omp critical
    {
      sum_approx_pi.num += indiv_approx_pi.num;
      sum_approx_pi.den += indiv_approx_pi.den;
    }
  }

  #pragma omp barrier

  #pragma omp single
  {
    double exp_pi_ratio = (double) sum_approx_pi.num / sum_approx_pi.den;
    printf("avg: %.10lf ~= %lu / %lu \n", exp_pi_ratio, sum_approx_pi.num, sum_approx_pi.den);
  }

  

  // best case
  /*

  |xxxx
  |    xx
  |     xx
  |       x
  |_______x___

  */

  /*
  #pragma omp single
  {
    Rational ideal_pi = calc_ideal_ratio();
  }
  */

  return 0;
}
