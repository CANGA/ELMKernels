
using ArrayI1 = Kokkos::View<int *>;
using ArrayB1 = Kokkos::View<bool *>;
using ArrayD1 = Kokkos::View<double *>;
using ArrayD2 = Kokkos::View<double **>;

struct CallCanFlux {

  CallCanFlux(ELM::LandType &Land_, ELM::VegProperties &Veg_, ArrayD2 &t_soisno_, ArrayD2 &h2osoi_ice_,
              ArrayD2 &h2osoi_liq_, ArrayD2 &dz_, ArrayD2 &rootfr_, ArrayD2 &sucsat_, ArrayD2 &watsat_, ArrayD2 &bsw_,
              ArrayD2 &smpso_, ArrayD2 &smpsc_, ArrayD2 &dleaf_, ArrayD2 &parsha_z_, ArrayD2 &parsun_z_,
              ArrayD2 &laisha_z_, ArrayD2 &laisun_z_, ArrayD2 &tlai_z_, ArrayD2 &rootr_, ArrayD2 &eff_porosity_,
              ArrayD2 &h2osoi_liqvol_, ArrayI1 &nrad_, ArrayI1 &snl_, ArrayI1 &frac_veg_nosno_, ArrayI1 &altmax_indx_,
              ArrayI1 &altmax_lastyear_indx_, ArrayD1 &frac_sno_, ArrayD1 &forc_hgt_u_patch_, ArrayD1 &thm_,
              ArrayD1 &thv_, ArrayD1 &max_dayl_, ArrayD1 &dayl_, ArrayD1 &tc_stress_, ArrayD1 &elai_, ArrayD1 &esai_,
              ArrayD1 &emv_, ArrayD1 &emg_, ArrayD1 &qg_, ArrayD1 &t_grnd_, ArrayD1 &forc_t_, ArrayD1 &forc_pbot_,
              ArrayD1 &forc_lwrad_, ArrayD1 &forc_u_, rrayD1 &forc_v_, ArrayD1 &forc_q_, ArrayD1 &forc_th_,
              ArrayD1 &z0mg_, ArrayD1 &dtime_, ArrayD1 &forc_hgt_t_patch_, ArrayD1 &forc_hgt_q_patch_, ArrayD1 &fwet_,
              ArrayD1 &fdry_, ArrayD1 &laisun_, ArrayD1 &laisha_, ArrayD1 &forc_rho_, ArrayD1 &snow_depth_,
              ArrayD1 &soilbeta_, ArrayD1 &frac_h2osfc_, ArrayD1 &t_h2osfc_, ArrayD1 &sabv_, ArrayD1 &h2ocan_,
              ArrayD1 &htop_, ArrayD1 &t10_, ArrayD1 &vcmaxcintsha_, ArrayD1 &vcmaxcintsun_, ArrayD1 &forc_pco2_,
              ArrayD1 &forc_po2_, ArrayD1 &qg_snow_, ArrayD1 &qg_soil_, ArrayD1 &qg_h2osfc_, ArrayD1 &dqgdT_,
              ArrayD1 &htvp_, ArrayD1 &eflx_sh_grnd_, ArrayD1 &eflx_sh_snow_, ArrayD1 &eflx_sh_soil_,
              ArrayD1 &eflx_sh_h2osfc_, ArrayD1 &qflx_evap_soi_, ArrayD1 &qflx_ev_snow_, ArrayD1 &qflx_ev_soil_,
              ArrayD1 &qflx_ev_h2osfc_, ArrayD1 &dlrad_, ArrayD1 &ulrad_, ArrayD1 &cgrnds_, ArrayD1 &cgrndl_,
              ArrayD1 &cgrnd_, ArrayD1 &t_ref2m_, ArrayD1 &t_ref2m_r_, ArrayD1 &q_ref2m_, ArrayD1 &rh_ref2m_,
              ArrayD1 &rh_ref2m_r_, ArrayD1 &qflx_tran_veg_, ArrayD1 &qflx_evap_veg_, ArrayD1 &eflx_sh_veg_,
              ArrayD1 &btran_, ArrayD1 &displa_, ArrayD1 &z0mv_, ArrayD1 &z0hv_, ArrayD1 &z0qv_, ArrayD1 &t_veg_)
      : Land(Land_), Veg(Veg_), t_soisno(t_soisno_), h2osoi_ice(h2osoi_ice_), h2osoi_liq(h2osoi_liq_), dz(dz_),
        rootfr(rootfr_), sucsat(sucsat_), watsat(watsat_), bsw(bsw_), smpso(smpso_), smpsc(smpsc_), dleaf(dleaf_),
        parsha_z(parsha_z_), parsun_z(parsun_z_), laisha_z(laisha_z_), laisun_z(laisun_z_), tlai_z(tlai_z_),
        rootr(rootr_), eff_porosity(eff_porosity_), h2osoi_liqvol(h2osoi_liqvol_), nrad(nrad_), snl(snl_),
        frac_veg_nosno(frac_veg_nosno_), altmax_indx(altmax_indx_), altmax_lastyear_indx(altmax_lastyear_indx_),
        frac_sno(frac_sno_), forc_hgt_u_patch(forc_hgt_u_patch_), thm(thm_), thv(thv_), max_dayl(max_dayl_),
        dayl(dayl_), tc_stress(tc_stress_), elai(elai_), esai(esai_), emv(emv_), emg(emg_), qg(qg_), t_grnd(t_grnd_),
        forc_t(forc_t_), forc_pbot(forc_pbot_), forc_lwrad(forc_lwrad_), forc_u(forc_u_), forc_v(forc_v_),
        forc_q(forc_q_), forc_th(forc_th_), z0mg(z0mg_), dtime(dtime_), forc_hgt_t_patch(forc_hgt_t_patch_),
        forc_hgt_q_patch(forc_hgt_q_patch_), fwet(fwet_), fdry(fdry_), laisun(laisun_), laisha(laisha_),
        forc_rho(forc_rho_), snow_depth(snow_depth_), soilbeta(soilbeta_), frac_h2osfc(frac_h2osfc_),
        t_h2osfc(t_h2osfc_), sabv(sabv_), h2ocan(h2ocan_), htop(htop_), t10(t10_), vcmaxcintsha(vcmaxcintsha_),
        vcmaxcintsun(vcmaxcintsun_), forc_pco2(forc_pco2_), forc_po2(forc_po2_), qg_snow(qg_snow_), qg_soil(qg_soil_),
        qg_h2osfc(qg_h2osfc_), dqgdT(dqgdT_), htvp(htvp_), eflx_sh_grnd(eflx_sh_grnd_), eflx_sh_snow(eflx_sh_snow_),
        eflx_sh_soil(eflx_sh_soil_), eflx_sh_h2osfc(eflx_sh_h2osfc_), qflx_evap_soi(qflx_evap_soi_),
        qflx_ev_snow(qflx_ev_snow_), qflx_ev_soil(qflx_ev_soil_), qflx_ev_h2osfc(qflx_ev_h2osfc_), dlrad(dlrad_),
        ulrad(ulrad_), cgrnds(cgrnds_), cgrndl(cgrndl_), cgrnd(cgrnd_), t_ref2m(t_ref2m_), t_ref2m_r(t_ref2m_r_),
        q_ref2m(q_ref2m_), rh_ref2m(rh_ref2m_), rh_ref2m_r(rh_ref2m_r_), qflx_tran_veg(qflx_tran_veg_),
        qflx_evap_veg(qflx_evap_veg_), eflx_sh_veg(eflx_sh_veg_), btran(btran_), displa(displa_), z0mv(z0mv_),
        z0hv(z0hv_), z0qv(z0qv_), t_veg(t_veg_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    double wtg;         // heat conductance for ground [m/s]
    double wtgq;        // latent heat conductance for ground [m/s]
    double wtalq;       // normalized latent heat cond. for air and leaf [-]
    double wtlq0;       // normalized latent heat conductance for leaf [-]
    double wtaq0;       // normalized latent heat conductance for air [-]
    double wtl0;        // normalized heat conductance for leaf [-]
    double wta0;        // normalized heat conductance for air [-]
    double wtal;        // normalized heat conductance for air and leaf [-]
    double dayl_factor; // scalar (0-1) for daylength effect on Vcmax
    double air;         // atmos. radiation temporay set
    double bir;         // atmos. radiation temporay set
    double cir;         // atmos. radiation temporay set
    double el;          // vapor pressure on leaf surface [pa]
    double qsatl;       // leaf specific humidity [kg/kg]
    double qsatldT;     // derivative of "qsatl" on "t_veg"
    double taf;         // air temperature within canopy space [K]
    double qaf;         // humidity of canopy air [kg/kg]
    double um;          // wind speed including the stablity effect [m/s]
    double ur;          // wind speed at reference height [m/s]
    double dth;         // diff of virtual temp. between ref. height and surface
    double dqh;         // diff of humidity between ref. height and surface
    double obu;         // Monin-Obukhov length (m)
    double zldis;       // reference height "minus" zero displacement height [m]
    double temp1;       // relation for potential temperature profile
    double temp2;       // relation for specific humidity profile
    double temp12m;     // relation for potential temperature profile applied at 2-m
    double temp22m;     // relation for specific humidity profile applied at 2-m
    double tlbef;       // leaf temperature from previous iteration [K]
    double delq;        // temporary
    double dt_veg;      // change in t_veg, last iteration (Kelvin)

    ELM::InitializeFlux_Can(
        const LandType &Land, const int &snl, const int &frac_veg_nosno, const double &frac_sno,
        const double &forc_hgt_u_patch, const double &thm, const double &thv, const double &max_dayl,
        const double &dayl, const int &altmax_indx, const int &altmax_lastyear_indx,
        const double t_soisno[nlevgrnd + nlevsno], const double h2osoi_ice[nlevgrnd + nlevsno],
        const double h2osoi_liq[nlevgrnd + nlevsno], const double dz[nlevgrnd + nlevsno], const double rootfr[nlevgrnd],
        const double &tc_stress, const double *sucsat, const double *watsat, const double *bsw, const double *smpso,
        const double *smpsc, const double &elai, const double &esai, const double &emv, const double &emg,
        const double &qg, const double &t_grnd, const double &forc_t, const double &forc_pbot, const double &forc_lwrad,
        const double &forc_u, const double &forc_v, const double &forc_q, const double &forc_th, const double &z0mg,

        double &btran, double &displa, double &z0mv, double &z0hv, double &z0qv,

        double *rootr, double eff_porosity[nlevgrnd], double h2osoi_liqvol[nlevgrnd + nlevsno], double &dayl_factor,
        double &air, double &bir, double &cir, double &el, double &qsatl, double &qsatldT, double &taf, double &qaf,
        double &um, double &ur, double &obu, double &zldis, double &delq, double &t_veg

    )
  }

private:
  ELM::LandType Land;
  ELM::VegProperties Veg;

  ArrayI1 nrad, snl, frac_veg_nosno, altmax_indx, altmax_lastyear_indx;

  ArrayD1 frac_sno, forc_hgt_u_patch, thm, thv, max_dayl, dayl, tc_stress, elai, esai, emv, emg, qg, t_grnd, forc_t,
      forc_pbot, forc_lwrad, forc_u, forc_v, forc_q, forc_th, z0mg, dtime, forc_hgt_t_patch, forc_hgt_q_patch, fwet,
      fdry, laisun, laisha, forc_rho, snow_depth, soilbeta, frac_h2osfc, t_h2osfc, sabv, h2ocan, htop, t10,
      vcmaxcintsha, vcmaxcintsun, forc_pco2, forc_po2, qg_snow, qg_soil, qg_h2osfc, dqgdT, htvp, eflx_sh_grnd,
      eflx_sh_snow, eflx_sh_soil, eflx_sh_h2osfc, qflx_evap_soi, qflx_ev_snow, qflx_ev_soil, qflx_ev_h2osfc, dlrad,
      ulrad, cgrnds, cgrndl, cgrnd, t_ref2m, t_ref2m_r, q_ref2m, rh_ref2m, rh_ref2m_r, qflx_tran_veg, qflx_evap_veg,
      eflx_sh_veg, btran, displa, z0mv, z0hv, z0qv, t_veg;

  ArrayD2 t_soisno, h2osoi_ice, h2osoi_liq, dz, rootfr, sucsat, watsat, bsw, smpso, smpsc, dleaf, parsha_z, parsun_z,
      laisha_z, laisun_z, tlai_z, rootr, eff_porosity, h2osoi_liqvol;
};

ELM::LandType &Land, ELM::VegProperties &Veg, ArrayD2 &t_soisno, ArrayD2 &h2osoi_ice, ArrayD2 &h2osoi_liq, ArrayD2 &dz,
    ArrayD2 &rootfr, ArrayD2 &sucsat, ArrayD2 &watsat, ArrayD2 &bsw, ArrayD2 &smpso, ArrayD2 &smpsc, ArrayD2 &dleaf,
    ArrayD2 &parsha_z, ArrayD2 &parsun_z, ArrayD2 &laisha_z, ArrayD2 &laisun_z, ArrayD2 &tlai_z, ArrayD2 &rootr,
    ArrayD2 &eff_porosity, ArrayD2 &h2osoi_liqvol, ArrayI1 &nrad, ArrayI1 &snl, ArrayI1 &frac_veg_nosno,
    ArrayI1 &altmax_indx, ArrayI1 &altmax_lastyear_indx, ArrayD1 &frac_sno, ArrayD1 &forc_hgt_u_patch, ArrayD1 &thm,
    ArrayD1 &thv, ArrayD1 &max_dayl, ArrayD1 &dayl, ArrayD1 &tc_stress, ArrayD1 &elai, ArrayD1 &esai, ArrayD1 &emv,
    ArrayD1 &emg, ArrayD1 &qg, ArrayD1 &t_grnd, ArrayD1 &forc_t, ArrayD1 &forc_pbot, ArrayD1 &forc_lwrad,
    ArrayD1 &forc_u, ArrayD1 &forc_v, ArrayD1 &forc_q, ArrayD1 &forc_th, ArrayD1 &z0mg, ArrayD1 &dtime,
    ArrayD1 &forc_hgt_t_patch, ArrayD1 &forc_hgt_q_patch, ArrayD1 &fwet, ArrayD1 &fdry, ArrayD1 &laisun,
    ArrayD1 &laisha, ArrayD1 &forc_rho, ArrayD1 &snow_depth, ArrayD1 &soilbeta, ArrayD1 &frac_h2osfc, ArrayD1 &t_h2osfc,
    ArrayD1 &sabv, ArrayD1 &h2ocan, ArrayD1 &htop, ArrayD1 &t10, ArrayD1 &vcmaxcintsha, ArrayD1 &vcmaxcintsun,
    ArrayD1 &forc_pco2, ArrayD1 &forc_po2, ArrayD1 &qg_snow, ArrayD1 &qg_soil, ArrayD1 &qg_h2osfc, ArrayD1 &dqgdT,
    ArrayD1 &htvp, ArrayD1 &eflx_sh_grnd, ArrayD1 &eflx_sh_snow, ArrayD1 &eflx_sh_soil, ArrayD1 &eflx_sh_h2osfc,
    ArrayD1 &qflx_evap_soi, ArrayD1 &qflx_ev_snow, ArrayD1 &qflx_ev_soil, ArrayD1 &qflx_ev_h2osfc, ArrayD1 &dlrad,
    ArrayD1 &ulrad, ArrayD1 &cgrnds, ArrayD1 &cgrndl, ArrayD1 &cgrnd, ArrayD1 &t_ref2m, ArrayD1 &t_ref2m_r,
    ArrayD1 &q_ref2m, ArrayD1 &rh_ref2m, ArrayD1 &rh_ref2m_r, ArrayD1 &qflx_tran_veg, ArrayD1 &qflx_evap_veg,
    ArrayD1 &eflx_sh_veg, ArrayD1 &btran, ArrayD1 &displa, ArrayD1 &z0mv, ArrayD1 &z0hv, ArrayD1 &z0qv, ArrayD1 &t_veg,

    Land_ Veg_ t_soisno_ h2osoi_ice_ h2osoi_liq_ dz_ rootfr_ sucsat_ watsat_ bsw_ smpso_ smpsc_ dleaf_ parsha_z_
        parsun_z_ laisha_z_ laisun_z_ tlai_z_ rootr_ eff_porosity_ h2osoi_liqvol_ nrad_ snl_ frac_veg_nosno_
            altmax_indx_ altmax_lastyear_indx_ frac_sno_ forc_hgt_u_patch_ thm_ thv_ max_dayl_ dayl_ tc_stress_ elai_
                esai_ emv_ emg_ qg_ t_grnd_ forc_t_ forc_pbot_ forc_lwrad_ forc_u_ forc_v_ forc_q_ forc_th_ z0mg_ dtime_
                    forc_hgt_t_patch_ forc_hgt_q_patch_ fwet_ fdry_ laisun_ laisha_ forc_rho_ snow_depth_ soilbeta_
                        frac_h2osfc_ t_h2osfc_ sabv_ h2ocan_ htop_ t10_ vcmaxcintsha_ vcmaxcintsun_ forc_pco2_ forc_po2_
                            qg_snow_ qg_soil_ qg_h2osfc_ dqgdT_ htvp_ eflx_sh_grnd_ eflx_sh_snow_ eflx_sh_soil_
                                eflx_sh_h2osfc_ qflx_evap_soi_ qflx_ev_snow_ qflx_ev_soil_ qflx_ev_h2osfc_ dlrad_ ulrad_
                                    cgrnds_ cgrndl_ cgrnd_ t_ref2m_ t_ref2m_r_ q_ref2m_ rh_ref2m_ rh_ref2m_r_
                                        qflx_tran_veg_ qflx_evap_veg_ eflx_sh_veg_ btran_ displa_ z0mv_ z0hv_ z0qv_
                                            t_veg_