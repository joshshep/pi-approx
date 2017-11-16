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
#include <random>
#include <iostream>
#include <string>
#include <vector>


typedef struct rational {
  uint64_t num;
  uint64_t den;
} Rational;




std::string nameForNumber (uint64_t number) {
  std::vector<std::string> ones {"","one", "two", "three", "four", "five", "six", "seven", "eight", "nine"};
  std::vector<std::string> teens{"ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen","sixteen", "seventeen", "eighteen", "nineteen"};
  std::vector<std::string> tens {"", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"};

  /*
  if (number < 10) {
    return ones[number];
  } else if (number < 20) {
    return teens [number - 10];
  } else if (number < 100) {
    return tens[number / 10] + ((number % 10 != 0) ? " " + nameForNumber(number % 10) : "");
  } else 
  */
  if (number < 1000) {
    return std::to_string(number); //nameForNumber(number / 100) + " hundred" + ((number % 100 != 0) ? " " + nameForNumber(number % 100) : "");
  } else if (number < 1000000) {
    return nameForNumber(number / 1000) + " thousand" + ((number % 1000 != 0) ? " " + nameForNumber(number % 1000) : "");
  } else if (number < 1000000000) {
    return nameForNumber(number / 1000000) + " million" + ((number % 1000000 != 0) ? " " + nameForNumber(number % 1000000) : "");
  } else if (number < 1000000000000) {
    return nameForNumber(number / 1000000000) + " billion" + ((number % 1000000000 != 0) ? " " + nameForNumber(number % 1000000000) : "");
  } else if (number < 1000000000000000) {
    return nameForNumber(number / 1000000000000) + " trillion" + ((number % 1000000000000 != 0) ? " " + nameForNumber(number % 1000000000000) : "");
  }
  return "error";
}



uint64_t bal_rand_r(uint32_t* seed) {
  int32_t randval = rand_r(seed);
  while (randval == (-1 << 31)) {
    randval = rand_r(seed);
  }
  return abs(randval);
}

uint64_t calc_ideal_ratio() {
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

  return rands_in_circ;
}

uint64_t calc_exper_ratio(uint32_t rank, uint64_t num_rands) {
  uint32_t seed = rank;
  uint64_t rands_in_circ = 0;
  const uint64_t circ_radius_sq = (((uint64_t) 1 << 31) - 1) * (((uint64_t) 1 << 31) - 1);
  uint64_t step = 100000000; // 100 million
  uint64_t i_iter = 0;
  uint32_t istep = 0;
  //printf("iterations = %lu ; step = %lu\n", num_rands, step);
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

  return rands_in_circ;
}

uint64_t calc_pi_float(int32_t rank, uint64_t num_rands) {
  int seed = rank;
  const double circ_radius = 1; // this value does *not* matter. As long as it isn't zero
  const double circ_radius_sq = circ_radius*circ_radius;
  std::mt19937_64 engine(seed);

  std::uniform_real_distribution<double> sqDist(-circ_radius, circ_radius);

  uint64_t rands_in_circ = 0;
  uint64_t step = 100000000; // 100 million
  uint64_t i_iter = 0;
  uint32_t istep = 0;
  while (i_iter < num_rands) {
    uint64_t step_stop = std::min(i_iter + step, num_rands);
    for (; i_iter < step_stop; ++i_iter) {
      double x_rand = sqDist(engine);
      double y_rand = sqDist(engine);
      double dist_2_origin = x_rand * x_rand + y_rand * y_rand;
      if (dist_2_origin <= circ_radius_sq) {
        ++rands_in_circ;
      }
    }
    #ifdef VERBOSE
    double exp_pi_ratio = (double) 4 * rands_in_circ / step_stop;

    printf("r%2u step %u: %.10lf ~= 4 * %lu / %lu \n", rank, istep, exp_pi_ratio, rands_in_circ, step_stop);
    #endif
    ++istep;
  }

  return rands_in_circ;
}

int main(int argc, char** argv, char** envp) {

  if (argc != 3) {
    printf("Usage: %s [# of iterations per thread] [# of threads]\n", argv[0]);
    return -1;
  }
  
  int num_threads = atoi(argv[2]);
  #pragma omp parallel num_threads(num_threads)
  assert( num_threads == omp_get_num_threads());


  srand(0);

  //uint32_t rank = 0;
  uint64_t total_num_rands = std::stoull(argv[1]);
  assert( total_num_rands % num_threads == 0 );
  uint64_t num_rands_per_thread = total_num_rands / num_threads;
  std::string num_name = nameForNumber(total_num_rands);

  #ifdef VERBOSE
  printf("# of threads: %d\n", num_threads);
  printf("total iterations: %lu (%s)\n", total_num_rands, num_name.c_str());
  #endif

  uint64_t num_sum = 0;
  #pragma omp parallel for reduction( + : num_sum )
  for (int tid=0; tid < num_threads; ++tid) {
    num_sum += calc_pi_float(tid, num_rands_per_thread);//calc_exper_ratio(tid, num_rands);
  }

  double exp_pi_ratio = (double) 4 * num_sum / total_num_rands;
  #ifdef VERBOSE
  printf("avg: %.10lf ~= %lu / %lu\n", exp_pi_ratio, 4 * num_sum, total_num_rands);
  #else
  printf("%.10lf\n", exp_pi_ratio);
  #endif

  

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
