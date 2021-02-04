// CallSurfaceRadiation.hh
#pragma once
#include "clm_constants.h"
#include "landtype.h"
#include "Kokkos_Core.hpp"
#include "SurfaceRadiation.h"

using ArrayD1 = Kokkos::View<double*>;
using ArrayI1 = Kokkos::View<int*>;
using ArrayD2 = Kokkos::View<double**>;


struct CallSurfRad {


  CallSurfRad(ELM::LandType& Land_, ArrayI1& nrad_, ArrayI1& snl_, ArrayD1& elai_, ArrayD1& snow_depth_, ArrayD1& fsr_, 
    ArrayD1& laisun_, ArrayD1& laisha_, ArrayD1& sabg_soil_, ArrayD1& sabg_snow_, ArrayD1& sabg_, ArrayD1& sabv_, 
    ArrayD1& fsa_, ArrayD2& tlai_z_, ArrayD2& fsun_z_, ArrayD2& forc_solad_, ArrayD2& forc_solai_, ArrayD2& fabd_sun_z_, 
    ArrayD2& fabd_sha_z_, ArrayD2& fabi_sun_z_, ArrayD2& fabi_sha_z_, ArrayD2& parsun_z_, ArrayD2& parsha_z_, 
    ArrayD2& laisun_z_, ArrayD2& laisha_z_, ArrayD2& sabg_lyr_, ArrayD2& ftdd_, ArrayD2& ftid_, ArrayD2& ftii_, ArrayD2& fabd_, 
    ArrayD2& fabi_, ArrayD2& albsod_, ArrayD2& albsoi_, ArrayD2& albsnd_hst_, ArrayD2& albsni_hst_, ArrayD2& albgrd_, 
    ArrayD2& albgri_, ArrayD2& flx_absdv_, ArrayD2& flx_absdn_, ArrayD2& flx_absiv_, ArrayD2& flx_absin_, ArrayD2& albd_, 
    ArrayD2& albi_) 
  : 
    Land(Land_), nrad(nrad_), snl(snl_), elai(elai_), snow_depth(snow_depth_), fsr(fsr_), laisun(laisun_), laisha(laisha_), 
    sabg_soil(sabg_soil_), sabg_snow(sabg_snow_), sabg(sabg_), sabv(sabv_), fsa(fsa_), tlai_z(tlai_z_), fsun_z(fsun_z_), 
    forc_solad(forc_solad_), forc_solai(forc_solai_), fabd_sun_z(fabd_sun_z_), fabd_sha_z(fabd_sha_z_), fabi_sun_z(fabi_sun_z_), 
    fabi_sha_z(fabi_sha_z_), parsun_z(parsun_z_), parsha_z(parsha_z_), laisun_z(laisun_z_), laisha_z(laisha_z_), sabg_lyr(sabg_lyr_), 
    ftdd(ftdd_), ftid(ftid_), ftii(ftii_), fabd(fabd_), fabi(fabi_), albsod(albsod_), albsoi(albsoi_), albsnd_hst(albsnd_hst_), 
    albsni_hst(albsni_hst_), albgrd(albgrd_), albgri(albgri_), flx_absdv(flx_absdv_), flx_absdn(flx_absdn_), flx_absiv(flx_absiv_), 
    flx_absin(flx_absin_), albd(albd_), albi(albi_) {}


  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    double trd[ELM::numrad]; // transmitted solar radiation: direct (W/m**2)
    double tri[ELM::numrad]; // transmitted solar radiation: diffuse (W/m**2)

    ELM::CanopySunShadeFractions(Land, nrad[i], elai[i], Kokkos::subview(tlai_z, i, Kokkos::ALL), Kokkos::subview(fsun_z, i, Kokkos::ALL), 
      Kokkos::subview(forc_solad, i, Kokkos::ALL), Kokkos::subview(forc_solai, i, Kokkos::ALL), Kokkos::subview(fabd_sun_z, i, Kokkos::ALL), 
      Kokkos::subview(fabd_sha_z, i, Kokkos::ALL), Kokkos::subview(fabi_sun_z, i, Kokkos::ALL), Kokkos::subview(fabi_sha_z, i, Kokkos::ALL), 
      Kokkos::subview(parsun_z, i, Kokkos::ALL), Kokkos::subview(parsha_z, i, Kokkos::ALL), Kokkos::subview(laisun_z, i, Kokkos::ALL), 
      Kokkos::subview(laisha_z, i, Kokkos::ALL), laisun[i], laisha[i]);

    ELM::SurfRadZeroFluxes(Land, sabg_soil[i], sabg_snow[i], sabg[i], sabv[i], fsa[i], Kokkos::subview(sabg_lyr, i, Kokkos::ALL));

    ELM::SurfRadAbsorbed(Land, snl[i], Kokkos::subview(ftdd, i, Kokkos::ALL), Kokkos::subview(ftid, i, Kokkos::ALL), 
      Kokkos::subview(ftii, i, Kokkos::ALL), Kokkos::subview(forc_solad, i, Kokkos::ALL), Kokkos::subview(forc_solai, i, Kokkos::ALL), 
      Kokkos::subview(fabd, i, Kokkos::ALL), Kokkos::subview(fabi, i, Kokkos::ALL), Kokkos::subview(albsod, i, Kokkos::ALL), 
      Kokkos::subview(albsoi, i, Kokkos::ALL), Kokkos::subview(albsnd_hst, i, Kokkos::ALL), Kokkos::subview(albsni_hst, i, Kokkos::ALL), 
      Kokkos::subview(albgrd, i, Kokkos::ALL), Kokkos::subview(albgri, i, Kokkos::ALL), sabv[i], fsa[i], sabg[i], sabg_soil[i], 
      sabg_snow[i], trd, tri);

    ELM::SurfRadLayers(Land, snl[i], sabg[i], sabg_snow[i], snow_depth[i], Kokkos::subview(flx_absdv, i, Kokkos::ALL), 
      Kokkos::subview(flx_absdn, i, Kokkos::ALL), Kokkos::subview(flx_absiv, i, Kokkos::ALL), Kokkos::subview(flx_absin, i, Kokkos::ALL), 
      trd, tri, Kokkos::subview(sabg_lyr, i, Kokkos::ALL));

    ELM::SurfRadReflected(Land, Kokkos::subview(albd, i, Kokkos::ALL), Kokkos::subview(albi, i, Kokkos::ALL), 
      Kokkos::subview(forc_solad, i, Kokkos::ALL), Kokkos::subview(forc_solai, i, Kokkos::ALL), fsr[i]);

  }

private:
  ELM::LandType Land;
  ArrayI1 nrad, snl;
  ArrayD1 elai, snow_depth, fsr, laisun, laisha, sabg_soil, sabg_snow, sabg, sabv, fsa;

  ArrayD2 tlai_z, fsun_z, forc_solad, forc_solai, fabd_sun_z, fabd_sha_z, fabi_sun_z, fabi_sha_z, parsun_z, parsha_z, 
  laisun_z, laisha_z, sabg_lyr, ftdd, ftid, ftii, fabd, fabi, albsod, albsoi, albsnd_hst, albsni_hst, albgrd, albgri, 
  flx_absdv, flx_absdn, flx_absiv, flx_absin, albd, albi;

};

void SurfaceRadiationInvoke(const int& ncells_, ELM::LandType& Land_, ArrayI1& nrad_, ArrayI1& snl_, ArrayD1& elai_, 
  ArrayD1& snow_depth_, ArrayD1& fsr_, ArrayD1& laisun_, ArrayD1& laisha_, ArrayD1& sabg_soil_, ArrayD1& sabg_snow_, 
  ArrayD1& sabg_, ArrayD1& sabv_, ArrayD1& fsa_, ArrayD2& tlai_z_, ArrayD2& fsun_z_, ArrayD2& forc_solad_, ArrayD2& forc_solai_, 
  ArrayD2& fabd_sun_z_, ArrayD2& fabd_sha_z_, ArrayD2& fabi_sun_z_, ArrayD2& fabi_sha_z_, ArrayD2& parsun_z_, ArrayD2& parsha_z_, 
  ArrayD2& laisun_z_, ArrayD2& laisha_z_, ArrayD2& sabg_lyr_, ArrayD2& ftdd_, ArrayD2& ftid_, ArrayD2& ftii_, ArrayD2& fabd_, 
  ArrayD2& fabi_, ArrayD2& albsod_, ArrayD2& albsoi_, ArrayD2& albsnd_hst_, ArrayD2& albsni_hst_, ArrayD2& albgrd_, 
  ArrayD2& albgri_, ArrayD2& flx_absdv_, ArrayD2& flx_absdn_, ArrayD2& flx_absiv_, ArrayD2& flx_absin_, ArrayD2& albd_, ArrayD2& albi_)
{

  CallSurfRad call_surfrad( Land_, nrad_, snl_, elai_, snow_depth_, fsr_, laisun_, laisha_, sabg_soil_, sabg_snow_, sabg_, sabv_, 
    fsa_, tlai_z_, fsun_z_, forc_solad_, forc_solai_, fabd_sun_z_, fabd_sha_z_, fabi_sun_z_, fabi_sha_z_, parsun_z_, parsha_z_, 
    laisun_z_, laisha_z_, sabg_lyr_, ftdd_, ftid_, ftii_, fabd_, fabi_, albsod_, albsoi_, albsnd_hst_, albsni_hst_, albgrd_, albgri_, 
    flx_absdv_, flx_absdn_, flx_absiv_, flx_absin_, albd_, albi_);

  Kokkos::parallel_for(ncells_, call_surfrad);
}
