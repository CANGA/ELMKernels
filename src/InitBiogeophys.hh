// clm_inst_biogeophys

#include "clm_constants.h"

void InitBiogeophys() {

  int glcmask = 0; // hardwired

  // snow water
  // Note: Glacier_mec columns are initialized with half the maximum snow cover.
  // This gives more realistic values of qflx_glcice sooner in the simulation
  // for columns with net ablation, at the cost of delaying ice formation
  // in columns with net accumulation.
  if (ltype == istice) {
    h2osno = h2osno_max;
  } else if (ltype == istice_mec || (ltype == istsoil && glcmask > 0)) {
    h2osno = 1.0 * h2osno_max;
  } else {
    h2osno_col = 0.0;
  }
  snow_depth_col = h2osno_col / bdsno;
}

/*
  ! Initialize urban constants

    call urbanparams_vars%Init(bounds_proc)

    ! Initialize ecophys constants

    call veg_vp%Init()

    ! Initialize soil order related constants

    call soilorderconInit()

    ! Initialize lake constants

    call LakeConInit()

    ! Initialize surface albedo constants

    call SurfaceAlbedoInitTimeConst(bounds_proc)

    ! Initialize vertical data components

    call initVertical(bounds_proc,               &
         snow_depth_col(begc:endc),              &
         urbanparams_vars%thick_wall(begl:endl), &
         urbanparams_vars%thick_roof(begl:endl))

    ! Initialize clm->drv and drv->clm data structures

    call atm2lnd_vars%Init( bounds_proc )
    call lnd2atm_vars%Init( bounds_proc )

    ! Initialize glc2lnd and lnd2glc even if running without create_glacier_mec_landunit,
    ! because at least some variables (such as the icemask) are referred to in code that
    ! is executed even when running without glc_mec.
    call glc2lnd_vars%Init( bounds_proc )
    call lnd2glc_vars%Init( bounds_proc )

    ! If single-column determine closest latitude and longitude

    if (single_column) then
       call getfil (fsurdat, locfn, 0)
       call shr_scam_getCloseLatLon(locfn, scmlat, scmlon, &
            closelat, closelon, closelatidx, closelonidx)
    end if

    ! Initialization of public data types

    call temperature_vars%init(bounds_proc,      &
         urbanparams_vars%em_roof(begl:endl),    &
         urbanparams_vars%em_wall(begl:endl),    &
         urbanparams_vars%em_improad(begl:endl), &
         urbanparams_vars%em_perroad(begl:endl))

    call canopystate_vars%init(bounds_proc)

    call soilstate_vars%init(bounds_proc)

    call waterstate_vars%init(bounds_proc,         &
         h2osno_col(begc:endc),                    &
         snow_depth_col(begc:endc),                &
         soilstate_vars%watsat_col(begc:endc, 1:), &
         temperature_vars%t_soisno_col(begc:endc, -nlevsno+1:) )


    call waterflux_vars%init(bounds_proc)

    call chemstate_vars%Init(bounds_proc)
    ! WJS (6-24-14): Without the following write statement, the assertion in
    ! energyflux_vars%init fails with pgi 13.9 on yellowstone. So for now, I'm leaving
    ! this write statement in place as a workaround for this problem.
    call energyflux_vars%init(bounds_proc, temperature_vars%t_grnd_col(begc:endc))

    call aerosol_vars%Init(bounds_proc)

    call frictionvel_vars%Init(bounds_proc)

    call lakestate_vars%Init(bounds_proc)

    call photosyns_vars%Init(bounds_proc)

    call soilhydrology_vars%Init(bounds_proc, nlfilename)

    call solarabs_vars%Init(bounds_proc)

    call surfalb_vars%Init(bounds_proc)

    call surfrad_vars%Init(bounds_proc)

    call dust_vars%Init(bounds_proc)

    call glc_diagnostics_vars%Init(bounds_proc)

    */
