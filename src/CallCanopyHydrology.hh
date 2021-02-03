// CallCanopyHydrology.hh
#include "clm_constants.h"
#include "landtype.h"
#include "Kokkos_Core.hpp"
#include "CanopyHydrology.h"

using ArrayD1 = Kokkos::View<double*>;
using ArrayI1 = Kokkos::View<int*>;
using ArrayD2 = Kokkos::View<double**>;

struct CallCanHydro {
  
  CallCanHydro(ELM::LandType& Land_, ArrayI1& frac_veg_nosno_, ArrayI1& n_irrig_steps_left_, ArrayI1& snl_, 
    double& dtime_, ArrayD1& forc_rain_, ArrayD1& forc_snow_, ArrayD1& elai_, ArrayD1& esai_, 
    ArrayD1& h2ocan_, ArrayD1& irrig_rate_, ArrayD1& qflx_irrig_, ArrayD1& 
    qflx_prec_grnd_, ArrayD1& qflx_snwcp_liq_, ArrayD1& qflx_snwcp_ice_, ArrayD1& qflx_snow_grnd_, 
    ArrayD1& qflx_rain_grnd_, ArrayD1& fwet_, ArrayD1& fdry_, ArrayD1& forc_t_, ArrayD1& t_grnd_, 
    ArrayD1& qflx_snow_melt_, ArrayD1& n_melt_, ArrayD1& micro_sigma_, ArrayD1& snow_depth_, ArrayD1& h2osno_, 
    ArrayD1& int_snow_, ArrayD1& qflx_snow_h2osfc_, ArrayD1& h2osfc_, ArrayD1& frac_h2osfc_, ArrayD1& frac_sno_eff_, ArrayD1& frac_sno_, ArrayD2& swe_old_, 
    ArrayD2& h2osoi_liq_, ArrayD2& h2osoi_ice_, ArrayD2& t_soisno_, ArrayD2& frac_iceold_, ArrayD2& dz_, ArrayD2& z_, 
    ArrayD2& zi_, ArrayD2& snw_rds_) 
  : 
    Land(Land_), frac_veg_nosno(frac_veg_nosno_), n_irrig_steps_left(n_irrig_steps_left_), snl(snl_), 
    dtime(dtime_), forc_rain(forc_rain_), forc_snow(forc_snow_), elai(elai_), esai(esai_), 
    h2ocan(h2ocan_), irrig_rate(irrig_rate_), qflx_irrig(qflx_irrig_), 
    qflx_prec_grnd(qflx_prec_grnd_), qflx_snwcp_liq(qflx_snwcp_liq_), qflx_snwcp_ice(qflx_snwcp_ice_), 
    qflx_snow_grnd(qflx_snow_grnd_), qflx_rain_grnd(qflx_rain_grnd_), fwet(fwet_), fdry(fdry_), forc_t(forc_t_), 
    t_grnd(t_grnd_), qflx_snow_melt(qflx_snow_melt_), n_melt(n_melt_), micro_sigma(micro_sigma_),
    snow_depth(snow_depth_), h2osno(h2osno_), int_snow(int_snow_), qflx_snow_h2osfc(qflx_snow_h2osfc_), h2osfc(h2osfc_), 
    frac_h2osfc(frac_h2osfc_), frac_sno_eff(frac_sno_eff_), frac_sno(frac_sno_), swe_old(swe_old_), h2osoi_liq(h2osoi_liq_), 
    h2osoi_ice(h2osoi_ice_), t_soisno(t_soisno_), frac_iceold(frac_iceold_), dz(dz_), z(z_), zi(zi_), snw_rds(snw_rds_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    static const double dewmx = 0.1; // hardwired
    static const bool do_capsnow = false; // hardwired
    static const int oldfflag = 1; // hardwired
    double qflx_candrip;
    double qflx_through_snow;
    double qflx_through_rain;
    double fracsnow;
    double fracrain;

    ELM::Interception(Land, frac_veg_nosno[i], forc_rain[i], forc_snow[i], dewmx, elai[i], esai[i], dtime, 
      h2ocan[i], qflx_candrip, qflx_through_snow, qflx_through_rain, fracsnow, fracrain);

    Irrigation(Land, irrig_rate[i], n_irrig_steps_left[i], qflx_irrig[i]);

    GroundFlux(Land, do_capsnow, frac_veg_nosno[i], forc_rain[i], forc_snow[i], qflx_irrig[i], qflx_candrip, 
      qflx_through_snow, qflx_through_rain, fracsnow, fracrain, qflx_prec_grnd[i], qflx_snwcp_liq[i], 
      qflx_snwcp_ice[i], qflx_snow_grnd[i], qflx_rain_grnd[i]);

    FracWet(Land, frac_veg_nosno[i], dewmx, elai[i], esai[i], h2ocan[i], fwet[i], fdry[i]);

    SnowInit(Land, dtime, do_capsnow, oldfflag, forc_t[i], t_grnd[i], qflx_snow_grnd[i], qflx_snow_melt[i], n_melt[i], snow_depth[i], 
      h2osno[i], int_snow[i], Kokkos::subview(swe_old, i, Kokkos::ALL), Kokkos::subview(h2osoi_liq, i, Kokkos::ALL),  
      Kokkos::subview(h2osoi_ice, i, Kokkos::ALL), Kokkos::subview(t_soisno, i, Kokkos::ALL), Kokkos::subview(frac_iceold, i, Kokkos::ALL), 
      snl[i], Kokkos::subview(dz, i, Kokkos::ALL), Kokkos::subview(z, i, Kokkos::ALL), Kokkos::subview(zi, i, Kokkos::ALL), 
      Kokkos::subview(snw_rds, i, Kokkos::ALL), qflx_snow_h2osfc[i], frac_sno_eff[i], frac_sno[i]);

    FracH2OSfc(Land, micro_sigma[i], h2osno[i], h2osfc[i], Kokkos::subview(h2osoi_liq, i, Kokkos::ALL), frac_sno[i], 
      frac_sno_eff[i], frac_h2osfc[i]);
  }

private:
  ELM::LandType Land;

  ArrayI1 frac_veg_nosno, n_irrig_steps_left, snl;

  ArrayD1 forc_rain, forc_snow, elai, esai, h2ocan, irrig_rate, qflx_irrig, 
  qflx_prec_grnd, qflx_snwcp_liq, qflx_snwcp_ice, qflx_snow_grnd, qflx_rain_grnd, fwet, fdry, forc_t, t_grnd, 
  qflx_snow_melt, n_melt, micro_sigma, snow_depth, h2osno, int_snow, qflx_snow_h2osfc, h2osfc, frac_h2osfc, frac_sno_eff, frac_sno;

  ArrayD2 swe_old, h2osoi_liq, h2osoi_ice, t_soisno, frac_iceold, dz, z, zi, snw_rds;

  double dtime;
};

void CanopyHydrologyInvoke(const int& ncells, ELM::LandType& Land_, ArrayI1& frac_veg_nosno_, ArrayI1& n_irrig_steps_left_, 
  ArrayI1& snl_, double& dtime_, ArrayD1& forc_rain_, ArrayD1& forc_snow_, ArrayD1& elai_, 
  ArrayD1& esai_, ArrayD1& h2ocan_, ArrayD1& irrig_rate_, 
  ArrayD1& qflx_irrig_, ArrayD1& qflx_prec_grnd_, ArrayD1& qflx_snwcp_liq_, ArrayD1& qflx_snwcp_ice_, 
  ArrayD1& qflx_snow_grnd_, ArrayD1& qflx_rain_grnd_, ArrayD1& fwet_, ArrayD1& fdry_, ArrayD1& forc_t_, 
  ArrayD1& t_grnd_, ArrayD1& qflx_snow_melt_, ArrayD1& n_melt_, ArrayD1& micro_sigma_, ArrayD1& snow_depth_, 
  ArrayD1& h2osno_, ArrayD1& int_snow_, ArrayD1& qflx_snow_h2osfc_, ArrayD1& h2osfc_, ArrayD1& frac_h2osfc_, 
  ArrayD1& frac_sno_eff_, ArrayD1& frac_sno_, ArrayD2& swe_old_, ArrayD2& h2osoi_liq_, ArrayD2& h2osoi_ice_, 
  ArrayD2& t_soisno_, ArrayD2& frac_iceold_, ArrayD2& dz_, ArrayD2& z_, ArrayD2& zi_, ArrayD2& snw_rds_) {

  CallCanHydro call_canhydro(Land_, frac_veg_nosno_, n_irrig_steps_left_, snl_, dtime_, forc_rain_, 
    forc_snow_, elai_, esai_, h2ocan_, irrig_rate_, qflx_irrig_, qflx_prec_grnd_, qflx_snwcp_liq_, 
    qflx_snwcp_ice_, qflx_snow_grnd_, qflx_rain_grnd_, fwet_, fdry_, forc_t_, t_grnd_, qflx_snow_melt_, 
    n_melt_, micro_sigma_, snow_depth_, h2osno_, int_snow_, qflx_snow_h2osfc_, h2osfc_, frac_h2osfc_, frac_sno_eff_, frac_sno_, 
    swe_old_, h2osoi_liq_, h2osoi_ice_, t_soisno_, frac_iceold_, dz_, z_, zi_, snw_rds_);

  Kokkos::parallel_for(ncells, call_canhydro);
}
