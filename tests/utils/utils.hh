//! A set of utilities for testing ELM kernels in C++

#ifndef ELM_KERNEL_TEST_UTILS_HH_
#define ELM_KERNEL_TEST_UTILS_HH_

#include <iostream>
#include <memory>
#include <type_traits>
#include <chrono>

#include "mpi.h"

#include "../utils/array.hh"

namespace ELM {
namespace Utils {

//
// Determine month_of_year from day_of_year
//
const static std::array<int,12> days_per_month = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
//
int month_from_day(int doy, int start_month) { //i_month will be the month inside the year, but we need the aboslute month from the starting month and year
  
  int doyM, total, i_month, ny;
  doyM = doy;
  for(int i=0;i<start_month-1;i++){
   doyM += days_per_month[i];
  }
  
  ny=doyM/365;
  //std::cout << "doyM: " << doyM << " ny: " << ny << std::endl;
  doyM=doyM%365;

  i_month = 0;
  total = 0;
  while (i_month < 12) {
    total += days_per_month[i_month];
    if (doyM < total) {
      break;
    } else {
      i_month++;
    }
  }
  assert(i_month < 12 && "Broken month_from_day.");
  return (ny*12 + i_month-start_month+1);
}


//
// Domain decomposition in x, y
//
std::pair<int, int>
get_domain_decomposition(int nprocs, int argc, char** argv) {
  	int numprocs_x=1;
	int numprocs_y=1;
	int n=nprocs;
	int i=2;
	while(n!=1){
		if(n%i==0)
		{
			if(numprocs_x-numprocs_y>=0){
				numprocs_y*=i;
			}else{
				numprocs_x*=i;
			}
			n=n/i;
		}else{
			i++;
		}
	}
  return std::make_pair(numprocs_x,numprocs_y);
}


//
// min/max/sum an array
//
template<size_t D>
std::array<double, 3>
min_max_sum(const MPI_Comm& comm, const Array<double,D>& arr)
{
  auto min_max = std::minmax_element(arr.begin(), arr.end());
  double min = *min_max.first;
  double max = *min_max.second;
  double sum = std::accumulate(arr.begin(), arr.end(), 0.);

  double gmin, gmax, gsum;
  MPI_Reduce(&min, &gmin, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  return std::array<double,3>{gmin, gmax, gsum};
}
  

// Convert precipitation to rain and snow
template<class Array_t>
void
convert_precip_to_rain_snow(Array_t& rain,
                            Array_t& snow, 
                            const Array_t& temp)
{
  const int nt = rain.extent(0);
  const int ng = rain.extent(1);

  for (int k=0; k<nt; k++) {
    for (int j=0; j<ng; j++) {
      if (temp(k,j) < 273.15) {
        // no need to update snow, both are initially total precip
        rain(k,j) = 0.0; 
      } else {
        // no need to update rain, both are initially total precip
        snow(k,j) = 0.0;
      }
    }
  }
}



namespace Clock {
using time_point_type = std::chrono::high_resolution_clock::time_point;
using duration_type = std::chrono::duration<double>;

time_point_type time() {
  return std::chrono::high_resolution_clock::now();
}

std::array<double,3> min_max_mean(const MPI_Comm& comm, duration_type duration) {
  double duration_d(duration.count());
  double min, max, mean;
  MPI_Reduce(&duration_d, &max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  MPI_Reduce(&duration_d, &min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&duration_d, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  mean /= numprocs;

  return std::array<double,3>{min, max, mean};
}

} // namespace Clock




} // namespace Utils
} // namespace ELM


#endif
