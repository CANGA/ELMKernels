//! A set of utilities for testing ELM kernels in C++
#ifndef ELM_KERNEL_TEST_NETCDF_HH_
#define ELM_KERNEL_TEST_NETCDF_HH_

#include <string>
#include <array>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "pnetcdf.h"


#define NCMPI_HANDLE_ERROR( status, what )         \
do {                                         \
   if ( status )                             \
   {                                         \
      std::cout                              \
         << __FILE__                         \
         << ':'                              \
         << __LINE__                         \
         << ':'                              \
         << what                             \
         << " failed with rc = "             \
         << status                           \
         << ':'                              \
         << ncmpi_strerror( status )            \
         << '\n' ;                           \
      abort() ;                              \
   }                                         \
} while ( 0 )





namespace ELM {
namespace Utils {

//
// Reads n_grid_cells worth of phenology data from NetCDF files, each with
// n_pft PFTs, into LAI and SAI matrices with a potential offset (in grid
// cells) given by offset.
// -----------------------------------------------------------------------------
template<typename Matrix_t>
void read_phenology(const std::string& fname,
                    int n_grid_cells, int n_pfts,
                    size_t offset,
                    Matrix_t& lai, Matrix_t& sai) {
	
  MPI_Info info;

  int ncid = -1;
  auto status = ncmpi_open(MPI_COMM_WORLD, fname.c_str(), NC_NOWRITE, info, &ncid);
  NCMPI_HANDLE_ERROR(status, std::string("nc_open")+" \""+fname+"\"");

  auto start = std::array<MPI_Offset,4>{0,0,0,0};
  auto count = std::array<MPI_Offset,4>{n_grid_cells, n_pfts, 1, 1};

  // get id and read LAI
  int varid = -1;
  status = ncmpi_inq_varid(ncid, "MONTHLY_LAI", &varid);
  NCMPI_HANDLE_ERROR(status, "nc_inq_varid for LAI");

  std::vector<double> data(n_grid_cells * n_pfts);
  status = ncmpi_get_vara_double_all(ncid, varid, start.data(), count.data(), data.data());
  NCMPI_HANDLE_ERROR(status, "nc_get_vara_double monthly_lai");

  // copy into array -- note this may thrash cache depending upon order of
  // LAI.  Good for C order, bad for Fortran
  for (int i=offset; i!=offset+n_grid_cells; ++i) {
    for (int j=0; j!=n_pfts; ++j) {
      lai[i][j] = data[(i-offset)*n_pfts + j];
    }
  }
  
  
  // get id and read SAI directly into the array
  status = ncmpi_inq_varid(ncid, "MONTHLY_SAI", &varid);
  NCMPI_HANDLE_ERROR( status, "nc_inq_varid for SAI" ) ;
  
  status = ncmpi_get_vara_double_all(ncid, varid, start.data(), count.data(), data.data());
  NCMPI_HANDLE_ERROR(status, "nc_get_vara_double monthly_sai");
  // copy into array -- note this may thrash cache depending upon order of
  // LAI.  Good for C order, bad for Fortran
  for (int i=offset; i!=offset+n_grid_cells; ++i) {
    for (int j=0; j!=n_pfts; ++j) {
      sai[i][j] = data[(i-offset)*n_pfts + j];
    }
  }

  status = ncmpi_close(ncid);
  NCMPI_HANDLE_ERROR( status, "nc_close" ) ;
}

void read_dimensions(const std::string& dir, std::size_t start_year, std::size_t start_month, std::size_t n_months){
	
	MPI_Info info;
	std::vector<int> all_lon(n_months*n_years); 
	std::vector<int> all_lat(n_months*n_years);

  	std::string basename="clmforc.GSWP3.c2011.0.5x0.5.Prec.";

	for (size_t yy=0; yy!=n_years; ++yy) {	
		for (size_t mm=0; mm!=n_months; ++mm) {
		    std::stringstream fname_full;
  			fname_full << dir << basename << start_year + yy << "-" << std::setw(2) << std::setfill('0') << mm+1 << ".nc";
			int indexunroll=yy*n_months+mm;

			int ncid = -1;
			int dimid = -1;
			//open
			auto status = ncmpi_open(MPI_COMM_WORLD,fname_full.str().c_str(), NC_NOWRITE, info, &ncid);
			NCMPI_HANDLE_ERROR(status, std::string("ncmpi_open")+" \""+fname_full.str().c_str()+"\"");

			//lon
			status = ncmpi_inq_dimid(ncid, "lon", &dimid);
			NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimid");
			MPI_Offset tmp_lon;
			status = ncmpi_inq_dimlen(ncid, dimid, &tmp_lon);
			NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimlen");
			all_lon[indexunroll]=tmp_lon;

			//lat
			status = ncmpi_inq_dimid(ncid, "lat", &dimid);
			NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimid");
			MPI_Offset tmp_lat;
			status = ncmpi_inq_dimlen(ncid, dimid, &tmp_lat);
			NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimlen");
			all_lat[indexunroll]=tmp_lat;

			//time
			status = ncmpi_inq_dimid(ncid, "time", &dimid);
			NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimid");
			MPI_Offset tmp_times;
			status = ncmpi_inq_dimlen(ncid, dimid, &tmp_times);
			NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimlen");
			n_total_times[indexunroll]=tmp_times;
			
			//close
			status = ncmpi_close(ncid);
			NCMPI_HANDLE_ERROR( status, "ncmpi_close" ) ;
		}

	}


	//check if all long and lat values are the same acroos the years and months
	for (size_t ii=1; ii!=n_months*n_years; ++ii){
		if(all_lon[ii] !=all_lon[ii-1] || all_lat[ii] !=all_lat[ii-1]){
			std::cout << "Error: All longitude or latitude values are not the same" << std::endl; 
			std::cout << "See year: " << ii/n_months << " month" << ii%n_months << std::endl; 
		}
	}

	//should be the same to the first component
	nx_glob=all_lon[0];
	ny_glob=all_lat[0];


}

//
// Reads NetCDF met forcing data, returning the min number of times read.
// -----------------------------------------------------------------------------
template<typename Matrix_t>
void read_forcing(const std::string& fname, std::size_t start_year, std::size_t start_months, std::size_t n_months, std::size_t i_beg, std::size_t j_beg) {
   
	MPI_Info info;

	size_t ntimes_l=0;
	for (size_t ii=0; ii!=n_months*n_years; ++ii) {
		ntimes_l+=n_total_times[ii];
	}
	size_t n_grid_cells=n_x_grid_cells*n_y_grid_cells;
	const int nf=2; //number of forcing data (in this case 2: precipitation and temperature)

	std::string basename[nf], varname[nf];
	std::vector<double> data_read(nf*ntimes_l*n_grid_cells);

	//precipitation
  	basename[0]="clmforc.GSWP3.c2011.0.5x0.5.Prec.";
  	varname[0]="PRECTmms";
	//temperature
	basename[1]="clmforc.GSWP3.c2011.0.5x0.5.TPQWL.";
	varname[1]="TBOT";

	int index_start=0;
	for (size_t yy=0; yy!=n_years; ++yy) {	
		for (size_t mm=0; mm!=n_months; ++mm) {
			int indexunroll=yy*n_months+mm;

			int ncid = -1;
			int dimid = -1;

			auto start = std::array<MPI_Offset,3>{0,j_beg,i_beg};
			auto count = std::array<MPI_Offset,3>{n_total_times[indexunroll], n_y_grid_cells, n_x_grid_cells};
			
			for(size_t nn=0;nn<nf;++nn){
				std::stringstream fname_full;
  				fname_full << dir << basename[nn] << start_year + yy << "-" << std::setw(2) << std::setfill('0') << mm+1 << ".nc";

				//open
				auto status = ncmpi_open(MPI_COMM_WORLD,fname_full.str().c_str(), NC_NOWRITE, info, &ncid);
				NCMPI_HANDLE_ERROR(status, std::string("ncmpi_open")+" \""+fname_full.str()+"\"");
				
				int varid = -1;
				
				//get id
				status = ncmpi_inq_varid(ncid,varname[nn].c_str(), &varid);
				NCMPI_HANDLE_ERROR(status, "nc_inq_varid");

				//read
				status = ncmpi_get_vara_double_all(ncid, varid, start.data(), count.data(), data_read.data()+index_start);
				NCMPI_HANDLE_ERROR( status, "nc_get_vara_double total_precip" );

				//close
				status = ncmpi_close(ncid);
				NCMPI_HANDLE_ERROR( status, "ncmpi_close" ) ;
				index_start+=n_total_times[indexunroll]*n_grid_cells;
			}
			
		}
	}
	
    for (int lcv_t=0; lcv_t<ntimes_l; ++lcv_t) {
	     for (int lcv_gc=0; lcv_gc<n_grid_cells; ++lcv_gc) {
		  		int indexunroll1=lcv_t*n_grid_cells+lcv_gc;
				//int indexunroll1=lcv_gc*ntimes_l+lcv_t; //not sure
				int indexunroll2=lcv_t*n_grid_cells+lcv_gc+n_grid_cells*ntimes_l;

      		if (data_read[indexunroll2] < 273.15) {
        			snow[lcv_t][lcv_gc] = data_read[indexunroll1];
        			rain[lcv_t][lcv_gc] = 0.;
        			temp[lcv_t][lcv_gc] = data_read[indexunroll2];
      		} else {
        			snow[lcv_t][lcv_gc] = 0.;
        			rain[lcv_t][lcv_gc] = data_read[indexunroll1];
        			temp[lcv_t][lcv_gc] = data_read[indexunroll2];
      		}
    	}
	}

}






} // namespace Utils
} // namespace ELM


#endif
