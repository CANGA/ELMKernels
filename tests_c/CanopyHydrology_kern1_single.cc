#include <netcdf.h>
#include <array>
#include <sstream>
#include <iterator>
#include <exception>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <iostream>
using namespace std;

#include "CanopyHydrology.hh"


#define handle_error( status, what )         \
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
         << nc_strerror( status )            \
         << '\n' ;                           \
      abort() ;                              \
   }                                         \
} while ( 0 )


int main(int argc, char ** argv)
{
	
	int ctype = 1;
	const int nmonths = 12;
	const int npfts = 17;
	bool urbpoi = false,do_capsnow=false;
	const double dtime = 1800.0;
	double elai = 0, esai = 0, dewmx = 0, forc_snow = 0, forc_rain = 0, qflx_prec_intr = 0, qflx_irrig = 0, qflx_prec_grnd = 0, qflx_snwcp_liq = 0, qflx_snwcp_ice = 0, qflx_snow_grnd_patch = 0, qflx_rain_grnd = 0,irrig_rate = 0 ,h2ocan = 0;
	int frac_veg_nosno = 1;
	int ipft,itime ;
	size_t ntimes=0;
	const int ltype = 1;
	double monthly_lai[nmonths][npfts] = {};
	double monthly_sai[nmonths][npfts] = {};
	int n_irrig_steps_left = 0;
	double * total_precip = NULL;
	
	int ncid= 0, varid= 0, dimid= 0, status= 0;
	

	/*! use netcdf*/
		
	/*! This will be the netCDF ID for the file and data variable.*/
	
	status = nc_open("links/surfacedataWBW.nc", NC_NOWRITE, &ncid);
	handle_error( status, "nc_create" ) ;
	
	
	status = nc_inq_varid(ncid, "MONTHLY_LAI", &varid);
	handle_error( status, "nc_inq_varid" ) ;


	array<size_t,4> start{ nmonths,npfts, 1, 1 };
	array<size_t,4> count{ 0,0,0,0 };
	
	status = nc_get_vara_double(ncid, varid, count.data(), start.data(),  *monthly_lai);
	handle_error( status, "nc_get_vara_double monthly_lai" ) ;
	
	
	
	
	status = nc_inq_varid(ncid, "MONTHLY_SAI", &varid);
	handle_error( status, "nc_inq_varid" ) ;
	/*!do ipft=1,npfts */
	/*! print *, monthly_lai(ipft,12), monthly_lai(ipft,3), monthly_lai(ipft,6) */
	/*!end do */
	ipft = 1;
	std::cout << "[ " << monthly_lai[12][ipft] << " , " << monthly_lai[3][ipft] << " , " << monthly_lai[6][ipft] << " ]" << std::endl;
	status = nc_get_vara_double(ncid, varid, count.data(), start.data(), *monthly_sai);
	handle_error( status, "nc_get_vara_double monthly_sai" ) ;
	
	status = nc_close(ncid);
	handle_error( status, "nc_close1" ) ;

	status = nc_open("links/forcing7.nc", NC_NOWRITE, &ncid);
	handle_error( status, "nc_create" ) ;

	status = nc_inq_dimid(ncid, "time", &dimid);
	handle_error( status, "nc_inq_dimid" ) ;
	

	status = nc_inq_dimlen(ncid, dimid, &ntimes );
	handle_error( status, "nc_inq_dimlen" ) ;

	array<size_t,3> start3{ntimes, 1, 1};
	array<size_t,3> count3{ 0,0,0 };
	total_precip = (double*)malloc(ntimes*1*1*sizeof(double));

	status = nc_inq_varid(ncid,"PRECTmms", &varid);
	handle_error( status, "nc_inq_varid" ) ;

	status = nc_get_vara_double(ncid, varid, count3.data(), start3.data(), total_precip);
	handle_error( status, "nc_get_vara_double total_precip" ) ;

	status = nc_close(ncid);
	handle_error( status, "nc_close2" ) ;

	//cout.precision(10);
	h2ocan = 0.0;
	for(itime = 0; itime < ntimes; itime += 1){
		forc_rain = total_precip[itime];
		forc_snow = 0.0;
		dewmx = 0.1;
		elai = monthly_lai[5][7];
		esai = monthly_sai[5][7];
		frac_veg_nosno = 1;
		irrig_rate = 0.0;
		n_irrig_steps_left = 0;
		urbpoi = false;
		do_capsnow = false;

                ELM::CanopyHydrologyKern1(dtime, forc_snow, forc_snow, irrig_rate,
                        ltype, ctype, urbpoi, do_capsnow,
                        elai, esai, dewmx, frac_veg_nosno,
                        h2ocan, n_irrig_steps_left,
                        qflx_prec_intr, qflx_irrig, qflx_prec_grnd,
                        qflx_snwcp_liq, qflx_snwcp_ice,
                        qflx_snow_grnd_patch, qflx_rain_grnd);
		
		//printf("[ %d,%E,%4.4f,%E,%E]\n",itime, forc_rain, h2ocan, qflx_prec_grnd, qflx_prec_intr);
		std::cout << "[ " << itime << " , " << forc_rain<< " , " << h2ocan<< " , " << qflx_prec_grnd<< " , " << qflx_prec_intr << " ]" << std::endl;
		
	}
	

	
	return 0;
}
