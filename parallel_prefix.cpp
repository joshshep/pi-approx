/*	Cpt S 483, Introduction to Parallel Computing
 *	School of Electrical Engineering and Computer Science
 *
 *  Joshua Shepherd
 *  2017-11-07
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

/*
Sync the processors (with barriers) and only print if rank==0
Otherwise acts as printf
*/
int printf_r0(int rank, const char* fmt, ...) {
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}

/*
Taken from lecture notes we have a matrix of the form
[ A 0 ]
[ B 1 ]

that we are taking the power of for the ith random number.

Here A and B are used repeatedly instead of the entire matrix
*/

/*
This is used in lieu of the full 2x2 matrix
*/
typedef struct pair {
  int A;
  int B;
} Pair;

class ParallelPrefix {
public:
  ParallelPrefix(int rank, int p, int series_len) {
    this->rank = rank;
    this->p = p;
    this->series_len = series_len;
    this->num_per_proc = series_len / p;
    if (rank == 0) { // we need an array to store the *whole* series
      buf = new int[series_len];
    } else {
      buf = new int[num_per_proc];
    }
    times = (rank == 0) ? new uint64_t[p] : NULL;
  
  }
  ~ParallelPrefix() {
    delete[] buf;
    delete[] times;
  }
  int basic_parallel_sum(int local_sum);
  int parallel_series();

  int simple_series();
  int simple_series2(Pair& mat);
  uint64_t time_simple_series();
  uint64_t time_parallel_series();
  uint64_t max_time_parallel_series();

  int matrix_iterate(Pair& mat, int n);
  int matrix_mult(Pair& mat);

  int print();
  int gather();
private:
  int rank, p;
  static const int a  = 17;
  static const int b  = 29;
  static const int P  = 47;
  static const int x0 = 5;
  int series_len;
  int num_per_proc;
  int* buf;
  uint64_t* times;
};

/*
Simple parallel prefix sum (not used)
*/
int ParallelPrefix::basic_parallel_sum(int local_sum) {
  int global_sum = local_sum;

  MPI_Status status;

  //printf("r%2d: starting sum = %2d\n", rank, local_sum);
  for (int order=1; order < this->p; order<<=1) {
    int neighbor_rank = this->rank ^ order;
    int neighbor_global_sum;

    //
    MPI_Sendrecv(&global_sum,          1, MPI_INT, neighbor_rank, 0, //send
                 &neighbor_global_sum, 1, MPI_INT, neighbor_rank, 0, //receive
                 MPI_COMM_WORLD, &status);
    global_sum += neighbor_global_sum;
    global_sum %= this->P;
    if (neighbor_rank < this->rank) {
      local_sum += neighbor_global_sum;
      local_sum %= this->P;
    }
  }
  //printf("r%2d: END local sum = %4d; global = %4d\n", rank, local_sum, global_sum);
  

  return local_sum;
}

/*
Performs parallel prefix on the matrices
*/
int ParallelPrefix::matrix_mult(Pair& mat) {
  Pair local = {1, 0}; // identity matrix
  Pair global = mat;

  MPI_Status status;

  //printf("r%2d: starting sum = %2d\n", rank, local_sum);
  for (int order=1; order < this->p; order<<=1) {
    int neighbor_rank = this->rank ^ order;
    Pair neighbor_global;

    //
    MPI_Sendrecv(&global,          sizeof(Pair), MPI_BYTE, neighbor_rank, 0, //send
                 &neighbor_global, sizeof(Pair), MPI_BYTE, neighbor_rank, 0, //receive
                 MPI_COMM_WORLD, &status);

    global.A = global.A * neighbor_global.A;
    global.B = global.B * neighbor_global.A + neighbor_global.B;
    global.A %= this->P;
    global.B %= this->P;
    if (neighbor_rank < this->rank) {
      local.A = local.A * neighbor_global.A; 
      local.B = local.B * neighbor_global.A + neighbor_global.B;
      local.A %= this->P;
      local.B %= this->P;
    }
  }

  mat.A = local.A;
  mat.B = local.B;
  //printf("r%2d: END local sum = %4d; global = %4d\n", rank, local_sum, global_sum);
  

  return 0;
}

int ParallelPrefix::matrix_iterate(Pair& mat, int n) {
  for (int i=0; i<n; ++i) {
    mat.A = mat.A * a;
    mat.B = a * mat.B + b;
    mat.A %= this->P;
    mat.B %= this->P;
  }
  return 0;
}  

/*
DEBUG this is used to gather the numbers from other processors in order to print them out sequentially
*/
int ParallelPrefix::gather() {
  MPI_Gather(buf, num_per_proc, MPI_INT, 
             buf, num_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
  return 0;
}

/*
DEBUG Print out the series of random numbers 
*/
int ParallelPrefix::print() {
  if (rank != 0) {
    return 0;
  }

  for (int i = 0; i < series_len; ++i) {
    printf("x_%8d: %5d\n", i, buf[i]);
  }
  return 0;
}

/*
Main processing
Compute the series of random numbers in parallel
*/
int ParallelPrefix::parallel_series() {
  Pair mat = {this->a, this->b};
  matrix_iterate(mat, num_per_proc - 1);
  matrix_mult(mat);
  simple_series2(mat);
  /*
  int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    MPI_Comm comm)*/

  return 0;
}  

/*
Serial test
Compute the series of random numbers in series
*/
int ParallelPrefix::simple_series() {
  int x = this->x0;
  int i = 0;
  buf[i] = x;
  for (i = 1; i < series_len; ++i) {
    x = (this->a * x + this->b) % this->P;
    buf[i] = x;
  }
  return 0;
}

/*
Same as simple_series() but can start from a Pair (i.e. not from the beginning)
*/
int ParallelPrefix::simple_series2(Pair& mat) {
  for (int i = 0; i < num_per_proc; ++i) {
    int x = (this->x0 * mat.A + mat.B) % this->P;
    buf[i] = x;
    matrix_iterate(mat, 1);
  }
  return 0;
}

// TODO how to pass a member function as a parameter?
// uint64_t ParallelPrefix::time_func( uint64_t (*f)()) {
uint64_t ParallelPrefix::time_simple_series() {
  struct timeval t1, t2;

  gettimeofday(&t1,NULL);
  simple_series();
  gettimeofday(&t2,NULL);
  
  uint64_t usecs = (t2.tv_sec-t1.tv_sec)*1e6 + (t2.tv_usec-t1.tv_usec);
  return usecs;
}

uint64_t ParallelPrefix::time_parallel_series() {
  struct timeval t1, t2;

  MPI_Barrier(MPI_COMM_WORLD);

  gettimeofday(&t1,NULL);
  parallel_series();
  gettimeofday(&t2,NULL);
  
  uint64_t usecs = (t2.tv_sec-t1.tv_sec)*1e6 + (t2.tv_usec-t1.tv_usec);
  return usecs;
}

/*
Find the maximum of the times of the parallel computation across all processes
*/
uint64_t ParallelPrefix::max_time_parallel_series() {
  uint64_t usecs = time_parallel_series();

  //gather from other processors
  MPI_Gather(&usecs, sizeof(uint64_t), MPI_BYTE, 
              times, sizeof(uint64_t), MPI_BYTE,
                0, MPI_COMM_WORLD);

  // find the max of the parallel times
  
  if (rank == 0) {
    int i_max = 0;
    for (int i=0; i < p; ++i) {
      if (times[i] > times[i_max])
        i_max = i;
    }
    return times[i_max];
  }
  return 0;
}


/*
Calculates the speedup

Note: this outputs in CSV for the report document
*/
int calc_speedup(int rank, int p, int series_len, int num_trials=100) {
  uint64_t usecs_serial_sum = 0;
  uint64_t usecs_parallel_sum = 0;  
  ParallelPrefix app(rank, p, series_len);
  
  //find serial time (omega)
  if (rank == 0) {
    for (int itrial=0; itrial<num_trials; ++itrial)
      usecs_serial_sum += app.time_simple_series();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  /*
  It is important to distinguish:
  We want the average of the max parallel times
  *NOT* the max of the average of the parallel times
  */
  // find parallel time T(n, p)
  for (int itrial=0; itrial<num_trials; ++itrial) {
    uint64_t usecs_max_parallel = app.max_time_parallel_series();

    // find the max of the parallel times
    usecs_parallel_sum += usecs_max_parallel;
  }

  if (rank != 0) {
    return 0;
  }

  // just r0

  // calculate speedup
  // note: num_trials cancels out as long as the sums are over the same number of trials
  double speedup = (double) usecs_serial_sum / usecs_parallel_sum;

  //output in csv-friendly format (`cat slurm*`)
  printf("%d,%lf\n", p, speedup);
  //printf("slowest: proc[%2d] = %lu usecs\n",i_max, times[i_max]);
  return 0;
}

/*
Calculates the runtimes for a variety of series lengths
Note: this outputs in CSV for the report document
*/
int gen_runtime_vs_n(int rank, int p, int num_trials=100) {
  // I care not about plural inconsistencies
  printf_r0(rank, "Series Length (n), %d Processes (p==%d)\n", p, p);
  for (int n=std::max(1,p); n <= 1<<20; n <<= 1) {
    uint64_t usecs_sum = 0;
    ParallelPrefix app(rank, p, n);
    for (int itrial=0; itrial<num_trials; ++itrial) {
      usecs_sum += app.max_time_parallel_series();
    }
    double usecs_avg = (double) usecs_sum / num_trials;
    printf_r0(rank, "%d,%lf\n", n, usecs_avg);

  }
  return 0;
}

/*
x_i = x_0 * a**i + b * (a**0 + a**1 + a**2 + ... + a**(i-1))
*/
int main(int argc, char** argv, char** envp) {

  int rank, p;
  int series_len = 1 << 20;
  //struct timeval t1, t2;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  /*
  ###########################################333
  Use calc_speedup() or gen_runtime_vs_n() but not both so you can `cat slurm*`
  ###########################################333
  */

  //speedup
  calc_speedup(rank, p, series_len);

  //runtime
  //gen_runtime_vs_n(rank, p);

  #if 0
  uint64_t usecs_serial, usecs_parallel;
  ParallelPrefix app(rank, p, series_len);
  if (rank == 0) {
    usecs_serial = app.time_simple_series();
    printf("N = %d\n", series_len);
    printf("Serial time: %lu usecs\n", usecs_serial);

    #ifdef DEBUG
    printf("~~~~~~~Serial~~~~~~~\n");
    app.print();
    printf("~~~~~~~~~~~~~~~~~~~~\n");
    #endif

  }
  // wait for r0 to catch up before we start in parallel
  MPI_Barrier(MPI_COMM_WORLD);


  usecs_parallel = app.time_parallel_series();

  #ifdef DEBUG
  app.gather();
  app.print();
  #endif

  //printf("r%2d: Parallel time: %lu usecs\n", rank, usecs_parallel);

  uint64_t* times = (rank == 0) ? new uint64_t[p] : NULL;
  
  MPI_Gather(&usecs_parallel, sizeof(uint64_t), MPI_BYTE, 
              times, sizeof(uint64_t), MPI_BYTE,
              0, MPI_COMM_WORLD);
  if (rank == 0) {
    uint64_t sum = 0;
    int i_max = 0;
    for (int i=0; i < p; ++i) {
      printf("p[%2d] = %lu usecs\n", i, times[i]);
      sum += times[i];
      if (times[i] > times[i_max])
        i_max = i;
    }
    double avg = (double) sum / p;
    printf("avg parallel time: %lf usecs\n", avg);
    printf("slowest: proc[%2d] = %lu usecs\n",i_max, times[i_max]);
    delete[] times;
  }
  /*
  printf("Finished?\n");
  //app.find_rands();
  */
  #endif 
  MPI_Finalize();
  return 0;
}
