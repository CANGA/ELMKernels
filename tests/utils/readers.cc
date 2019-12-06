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
#include "array.hh"


#define NCMPI_HANDLE_ERROR( status, what )      \
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
          << ncmpi_strerror( status )           \
          << '\n' ;                             \
      abort() ;                                 \
    }                                           \
  } while ( 0 )





namespace ELM {
namespace IO {

//
// Reads nx*ny worth of phenology data from NetCDF files, each with
// n_pft PFTs, into LAI and SAI matrices // -----------------------------------------------------------------------------
void read_phenology(const MPI_Comm& comm,
                    const std::string& dir,const std::string& basename, const std::string& phenology_type,
                    int start_year, int start_month,
                    int i_beg, int j_beg, ELM::Utils::Array<double,4>& arr){

  MPI_Info info;
  std::string varname;

	MPI_Info_create (&info);

  if(phenology_type.compare("ELAI")){
    varname="MONTHLY_LAI";
  }else{
    if(phenology_type.compare("ESAI")){
      varname="MONTHLY_SAI";
    }else{
      std::cout << "Error: phenology_type "<< phenology_type.c_str() << "doesn't exist" << std::endl; 
    }	
  }

  const auto dims=arr.shape();
  const int n_months = std::get<0>(dims);
  const int nx = std::get<1>(dims);
  const int ny = std::get<2>(dims);
  const int n_pfts = std::get<0>(dims);

  std::vector<double> data_read(n_months*nx*ny*n_pfts);

  int index_start=0; //index to track the position in the vector data_read
  for (int mm=0; mm!=n_months; ++mm) {
    int month=(start_month+mm-1)%12+1;
    int year=start_year+(start_month+mm-1)/12;
		
    //change this
    std::stringstream fname_full;
    fname_full << dir << basename << year << "-" << std::setw(2) << std::setfill('0') << month << ".nc";

    int ncid = -1;
    auto status = ncmpi_open(MPI_COMM_WORLD, fname_full.str().c_str(), NC_NOWRITE, info, &ncid);
    NCMPI_HANDLE_ERROR(status, std::string("nc_open")+" \""+fname_full.str().c_str()+"\"");

    auto start = std::array<MPI_Offset,4>{0,0,j_beg,i_beg};
    auto count = std::array<MPI_Offset,4>{1,n_pfts, ny, nx};

    // get id and read LAI
    int varid = -1;
    status = ncmpi_inq_varid(ncid, "MONTHLY_LAI", &varid);
    NCMPI_HANDLE_ERROR(status, "nc_inq_varid for LAI");

    status = ncmpi_get_vara_double_all(ncid, varid, start.data(), count.data(), data_read.data());
    NCMPI_HANDLE_ERROR(status, "nc_get_vara_double monthly_lai");

    status = ncmpi_close(ncid);
    NCMPI_HANDLE_ERROR( status, "nc_close" ) ;

    index_start+=n_pfts*nx*ny;
  }

  for(int ii=0;ii<n_months*n_pfts*nx*ny;ii++){
    //unpacking the indexes from 1D (ii) to 3D (k,l,i,j)
    int i,j,k,l; //k for time, l for pfts
    k = (ii/(nx*ny*n_pfts))     ;
    l = (ii/(nx*ny)) % n_pfts;
    j = (ii/(nx)) % ny;
    i = (ii) % nx;
    arr[k][l][i][j]=data_read[ii];
  }

}

std::tuple<int, int, int> 
get_dimensions(const std::string& dir, const std::string& basename, 
               int start_year, int start_month, int n_months){
	
  MPI_Info info;
  std::vector<int> all_lon(n_months); 
  std::vector<int> all_lat(n_months);

	MPI_Info_create (&info);


  int ntimes=0;
  for (int mm=0; mm!=n_months; ++mm) {
    int month=(start_month+mm-1)%12+1;
    int year=start_year+(start_month+mm-1)/12;

    std::stringstream fname_full;
    //the format should be like that
    fname_full << dir << basename << year << "-" << std::setw(2) << std::setfill('0') << month << ".nc";

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
    all_lon[mm]=tmp_lon;

    //lat
    status = ncmpi_inq_dimid(ncid, "lat", &dimid);
    NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimid");
    MPI_Offset tmp_lat;
    status = ncmpi_inq_dimlen(ncid, dimid, &tmp_lat);
    NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimlen");
    all_lat[mm]=tmp_lat;

    //time
    status = ncmpi_inq_dimid(ncid, "time", &dimid);
    NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimid");
    MPI_Offset tmp_times;
    status = ncmpi_inq_dimlen(ncid, dimid, &tmp_times);
    NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimlen");
    ntimes+=tmp_times;
			
    //close
    status = ncmpi_close(ncid);
    NCMPI_HANDLE_ERROR( status, "ncmpi_close" ) ;

  }

  //check if all long and lat values are the same acroos the years and months
  for (int mm=1; mm!=n_months; ++mm){
    if(all_lon[mm] !=all_lon[mm-1] || all_lat[mm] !=all_lat[mm-1]){
      std::cout << "Error: All longitude or latitude values are not the same" << std::endl; 
      std::cout << "See year: " << mm/n_months << " month" << mm%n_months << std::endl; 
    }
  }


//should be the same to the first component
  return std::make_tuple(ntimes, all_lon[0], all_lat[0]);

}

//
// Reads NetCDF met forcing data, returning the min number of times read.
// -----------------------------------------------------------------------------
void read_forcing(const MPI_Comm& comm,
                  const std::string& dir,const std::string& basename, const std::string& forcing_type,
                  int start_year, int start_month, int n_months,
                  int i_beg, int j_beg, ELM::Utils::Array<double,3>& arr){


  MPI_Info info;
	MPI_Info_create (&info);

  const auto dims=arr.shape();
  const int n_times = std::get<0>(dims);
  const int nx = std::get<1>(dims);
  const int ny = std::get<2>(dims);

  std::string varname;
	
  if(forcing_type.compare("PRECIP")){
    varname="PRECTmms";
  }else{
    if(forcing_type.compare("AIR_TEMP")){
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
    fname_full << dir << basename << year << "-" << std::setw(2) << std::setfill('0') << month << ".nc";

    int ncid = -1;
    int dimid = -1;

    //open
    auto status = ncmpi_open(MPI_COMM_WORLD,fname_full.str().c_str(), NC_NOWRITE, info, &ncid);
    NCMPI_HANDLE_ERROR(status, std::string("ncmpi_open")+" \""+fname_full.str()+"\"");
		
    //time
    status = ncmpi_inq_dimid(ncid, "time", &dimid);
    NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimid");
    MPI_Offset tmp_times;
    status = ncmpi_inq_dimlen(ncid, dimid, &tmp_times);
    NCMPI_HANDLE_ERROR(status, "ncmpi_inq_dimlen");

    int varid = -1;

    //define offsets
    auto start = std::array<MPI_Offset,3>{0,j_beg,i_beg};
    auto count = std::array<MPI_Offset,3>{tmp_times, ny, nx};

    //get id
    status = ncmpi_inq_varid(ncid,varname.c_str(), &varid);
    NCMPI_HANDLE_ERROR(status, "nc_inq_varid");

    //read
    status = ncmpi_get_vara_double_all(ncid, varid, start.data(), count.data(), (double*) arr[index_start].begin());
    NCMPI_HANDLE_ERROR( status, "nc_get_vara_double total_precip" );

    //close
    status = ncmpi_close(ncid);
    NCMPI_HANDLE_ERROR( status, "ncmpi_close" ) ;

    index_start+=tmp_times;
			
  }
}


void convert_precip_to_rain_snow(ELM::Utils::Array<double,3>& rain, ELM::Utils::Array<double,3>& snow, 
        ELM::Utils::Array<double,3>& temp){

  const auto dims=rain.shape();
  const int nt = std::get<0>(dims);
  const int nx = std::get<1>(dims);
  const int ny = std::get<2>(dims);

  for(int k=0;k<nt;k++){
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
        if(temp[k][i][j]<273.15){
          rain[k][i][j]=0.0; //no need to update the snow (it will have the correct value)
        }else{
          snow[k][i][j]=0.0; //no need to update the rain (it will have the correct value)
        }
      }
    }
  }

}

} // namespace Utils
} // namespace ELM


#endif
