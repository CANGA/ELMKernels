
#include "invoke_kernel.hh"
#include "atm_data.h"
#include "surface_radiation.h"
#include "surface_radiation_kokkos.hh"

void ELM::kokkos_surface_radiation(ELMStateType& S)
{
  const int ncols = S.num_columns;

  // local work arrays
  ViewD2 trd("trd", ncols, ELMdims::numrad());
  ViewD2 tri("tri", ncols, ELMdims::numrad());
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  // call surface_radiation kernels
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  /* =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
  auto surfrad_kernels = ELM_LAMBDA (const int& idx) {

    // call canopy_sunshade_fractions kernel
    ELM::surface_radiation::canopy_sunshade_fractions(
        S.Land,
        S.nrad(idx),
        S.elai(idx),
        Kokkos::subview(S.tlai_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fsun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.forc_solad, idx, Kokkos::ALL),
        Kokkos::subview(S.forc_solai, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd_sha_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi_sha_z, idx, Kokkos::ALL),
        Kokkos::subview(S.parsun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.parsha_z, idx, Kokkos::ALL),
        Kokkos::subview(S.laisun_z, idx, Kokkos::ALL),
        Kokkos::subview(S.laisha_z, idx, Kokkos::ALL),
        S.laisun(idx),
        S.laisha(idx));

    ELM::surface_radiation::initialize_flux(
        S.Land,
        S.sabg_soil(idx),
        S.sabg_snow(idx),
        S.sabg(idx),
        S.sabv(idx),
        S.fsa(idx),
        Kokkos::subview(S.sabg_lyr, idx, Kokkos::ALL));

    ELM::surface_radiation::total_absorbed_radiation(
        S.Land,
        S.snl(idx),
        Kokkos::subview(S.ftdd, idx, Kokkos::ALL),
        Kokkos::subview(S.ftid, idx, Kokkos::ALL),
        Kokkos::subview(S.ftii, idx, Kokkos::ALL),
        Kokkos::subview(S.forc_solad, idx, Kokkos::ALL),
        Kokkos::subview(S.forc_solai, idx, Kokkos::ALL),
        Kokkos::subview(S.fabd, idx, Kokkos::ALL),
        Kokkos::subview(S.fabi, idx, Kokkos::ALL),
        Kokkos::subview(S.albsod, idx, Kokkos::ALL),
        Kokkos::subview(S.albsoi, idx, Kokkos::ALL),
        Kokkos::subview(S.albsnd, idx, Kokkos::ALL),
        Kokkos::subview(S.albsni, idx, Kokkos::ALL),
        Kokkos::subview(S.albgrd, idx, Kokkos::ALL),
        Kokkos::subview(S.albgri, idx, Kokkos::ALL),
        S.sabv(idx),
        S.fsa(idx),
        S.sabg(idx),
        S.sabg_soil(idx),
        S.sabg_snow(idx),
        Kokkos::subview(trd, idx, Kokkos::ALL),
        Kokkos::subview(tri, idx, Kokkos::ALL));

    ELM::surface_radiation::layer_absorbed_radiation(
        S.Land,
        S.snl(idx),
        S.sabg(idx),
        S.sabg_snow(idx),
        S.snow_depth(idx),
        Kokkos::subview(S.flx_absdv, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absdn, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absiv, idx, Kokkos::ALL),
        Kokkos::subview(S.flx_absin, idx, Kokkos::ALL),
        Kokkos::subview(trd, idx, Kokkos::ALL),
        Kokkos::subview(tri, idx, Kokkos::ALL),
        Kokkos::subview(S.sabg_lyr, idx, Kokkos::ALL));

    ELM::surface_radiation::reflected_radiation(
        S.Land,
        Kokkos::subview(S.albd, idx, Kokkos::ALL),
        Kokkos::subview(S.albi, idx, Kokkos::ALL),
        Kokkos::subview(S.forc_solad, idx, Kokkos::ALL),
        Kokkos::subview(S.forc_solai, idx, Kokkos::ALL),
        S.fsr(idx));
  }; // end surfrad lambda
  apply_parallel_for(surfrad_kernels, "kokkos_surface_radiation", ncols);
}

