//! A set of utilities for testing ELM kernels in C++
#ifndef ELM_KERNEL_TEST_NETCDF_HH_
#define ELM_KERNEL_TEST_NETCDF_HH_

#include <string>
#include <array>
#include <vector>
#include <sstream>
#include <iostream>
#include <tuple>
#include <iomanip>
#include "netcdf.h"
#include "array.hh"


#define NC_HANDLE_ERROR( status, what )      \
  do {                                          \
    if ( status )                               \
    {                                           \
      std::cout                                 \
          << __FILE__                           \
          << ':'                                \
          << __LINE__                           \
          << ':'                                \
          << what                               \
          << " failed with rc = "               \
          << status                             \
          << ':'                                \
          << nc_strerror( status )           \
          << '\n' ;                             \
      abort() ;                                 \
    }                                           \
  } while ( 0 )





namespace ELM {
namespace IO {

//
// Reads nx*ny worth of phenology data from NetCDF files, each with
// n_pft PFTs, into LAI and SAI matrices // -----------------------------------------------------------------------------
void read_phenology(const std::string& dir,const std::string& basename, const std::string& phenology_type,
                    int start_year, int start_month,
                    ELM::Utils::Array<double,4>& arr){

  std::string varname;

  if(!phenology_type.compare("ELAI")){
    varname="MONTHLY_LAI";
  }else{
    if(!phenology_type.compare("ESAI")){
      varname="MONTHLY_SAI";
    }else{
      std::cout << "Error: phenology_type "<< phenology_type.c_str() << "doesn't exist" << std::endl; 
    }	
  }

  const auto dims=arr.shape();
  const int n_months = std::get<0>(dims);
  const int n_pfts = std::get<1>(dims);
  const int ny = std::get<2>(dims);
  const int nx = std::get<3>(dims);
  std::stringstream fname_full;
  fname_full << dir << "/" << basename;

  for (int mm=0; mm!=n_months; ++mm) {	
    int ncid = -1;
    auto status = nc_open(fname_full.str().c_str(), NC_NOWRITE, &ncid);
    NC_HANDLE_ERROR(status, std::string("nc_open")+" \""+fname_full.str().c_str()+"\"");

    auto start = std::array<size_t,4>{0,0,0,0};
    auto count = std::array<size_t,4>{1,n_pfts, ny, nx};

    // get id and read varname
    int varid = -1;
    status = nc_inq_varid(ncid, varname.c_str(), &varid);
	 NC_HANDLE_ERROR(status, "nc_inq_varid phenology");

    status = nc_get_vara_double(ncid, varid, start.data(), count.data(), (double*) arr[mm].begin());
    NC_HANDLE_ERROR(status, "nc_get_vara_double phenology");

    status = nc_close(ncid);
    NC_HANDLE_ERROR( status, "nc_close" ) ;
  }
}

std::tuple<int, int, int> 
get_dimensions(const std::string& dir, const std::string& basename, 
               int start_year, int start_month, int n_months){
	
  std::vector<int> all_lon(n_months); 
  std::vector<int> all_lat(n_months);


  int ntimes=0;
  for (int mm=0; mm!=n_months; ++mm) {
    int month=(start_month+mm-1)%12+1;
    int year=start_year+(start_month+mm-1)/12;

    std::stringstream fname_full;
    //the format should be like that
    fname_full << dir << "/" << basename << year << "-" << std::setw(2) << std::setfill('0') << month << ".nc";

    int ncid = -1;
    int dimid = -1;
    //open
    auto status = nc_open(fname_full.str().c_str(), NC_NOWRITE, &ncid);
    NC_HANDLE_ERROR(status, std::string("nc_open")+" \""+fname_full.str().c_str()+"\"");

    //lon
    status = nc_inq_dimid(ncid, "lon", &dimid);
    NC_HANDLE_ERROR(status, "nc_inq_dimid");
    size_t tmp_lon;
    status = nc_inq_dimlen(ncid, dimid, &tmp_lon);
    NC_HANDLE_ERROR(status, "nc_inq_dimlen");
    all_lon[mm]=tmp_lon;

    //lat
    status = nc_inq_dimid(ncid, "lat", &dimid);
    NC_HANDLE_ERROR(status, "nc_inq_dimid");
    size_t tmp_lat;
    status = nc_inq_dimlen(ncid, dimid, &tmp_lat);
    NC_HANDLE_ERROR(status, "nc_inq_dimlen");
    all_lat[mm]=tmp_lat;

    //time
    status = nc_inq_dimid(ncid, "time", &dimid);
    NC_HANDLE_ERROR(status, "nc_inq_dimid");
    size_t tmp_times;
    status = nc_inq_dimlen(ncid, dimid, &tmp_times);
    NC_HANDLE_ERROR(status, "nc_inq_dimlen");
    ntimes+=tmp_times;
			
    //close
    status = nc_close(ncid);
    NC_HANDLE_ERROR( status, "nc_close" ) ;

  }

  //check if all long and lat values are the same acroos the years and months
  for (int mm=1; mm!=n_months; ++mm){
    if(all_lon[mm] !=all_lon[mm-1] || all_lat[mm] !=all_lat[mm-1]){
      std::cout << "Error: All longitude or latitude values are not the same" << std::endl; 
      std::cout << "See year: " << mm/n_months << " month" << mm%n_months << std::endl; 
    }
  }


//should be the same to the first component
  return std::make_tuple(ntimes, all_lat[0], all_lon[0]);

}

//
// Reads NetCDF met forcing data, returning the min number of times read.
// -----------------------------------------------------------------------------
void read_forcing(const std::string& dir,const std::string& basename, const std::string& forcing_type,
                  int start_year, int start_month, int n_months,
                  ELM::Utils::Array<double,3>& arr){


  const auto dims=arr.shape();
  const int ny = std::get<1>(dims);
  const int nx = std::get<2>(dims);


  std::string varname;
	
  if(!forcing_type.compare("PRECIP")){
    varname="PRECTmms";
  }else{
    if(!forcing_type.compare("AIR_TEMP")){
      varname="TBOT";
    }else{
      std::cout << "Error: forcing_type "<< forcing_type.c_str() << "doesn't exist" << std::endl; 
    }	
  }

  int index_start=0; //index to track the position in the vector data_read
  for (int mm=0; mm!=n_months; ++mm) {
    int month=(start_month+mm-1)%12+1;
    int year=start_year+(start_month+mm-1)/12;

    std::stringstream fname_full;
    fname_full << dir << "/" << basename << year << "-" << std::setw(2) << std::setfill('0') << month << ".nc";

    int ncid = -1;
    int dimid = -1;

    //open
    auto status = nc_open(fname_full.str().c_str(), NC_NOWRITE, &ncid);
    NC_HANDLE_ERROR(status, std::string("nc_open")+" \""+fname_full.str()+"\"");
		
    //time
    status = nc_inq_dimid(ncid, "time", &dimid);
    NC_HANDLE_ERROR(status, "nc_inq_dimid");
	 size_t tmp_times;
    status = nc_inq_dimlen(ncid, dimid, &tmp_times);
    NC_HANDLE_ERROR(status, "nc_inq_dimlen");

    int varid = -1;

    //define offsets
    auto start = std::array<size_t,3>{0,0,0};
    auto count = std::array<size_t,3>{tmp_times, ny, nx};

    //get id
    status = nc_inq_varid(ncid,varname.c_str(), &varid);
    NC_HANDLE_ERROR(status, "nc_inq_varid forcing");

    //read
    status = nc_get_vara_double(ncid, varid, start.data(), count.data(), (double*) arr[index_start].begin());
    NC_HANDLE_ERROR( status, "nc_get_vara_double forcing" );

    //close
    status = nc_close(ncid);
    NC_HANDLE_ERROR( status, "nc_close" ) ;

    index_start+=tmp_times;
			
  }
}


void convert_precip_to_rain_snow(ELM::Utils::Array<double,2>& rain, ELM::Utils::Array<double,2>& snow, 
        ELM::Utils::Array<double,2>& temp){

  const auto dims=rain.shape();
  const int nt = std::get<0>(dims);
  const int ny = std::get<1>(dims);
  const int nx = std::get<2>(dims);


  for(int k=0;k<nt;k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
        if(temp[k][j][i]<273.15){
          rain[k][j][i]=0.0; //no need to update snow (it will have the correct value)
        }else{
          snow[k][j][i]=0.0; //no need to update rain (it will have the correct value)
        }
      }
    }
  }

}

} // namespace Utils
} // namespace ELM


#endif
