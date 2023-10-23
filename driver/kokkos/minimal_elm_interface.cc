
#include "minimal_elm_interface.hh"

using namespace ELM::ELMdims;


void ELM::MinimalInterface::setup(int ncols)
{ 
  Kokkos::initialize();
  ncols_ = ncols;
  S_ = std::make_shared<ELMStatType>(ncols_);
}

void ELM::MinimalInterface::initialize()
{
  /* need data from
  
  constants:
    grid:
      dz
      zsoi - centroids
      zisoi - interfaces

    soil:
      isoicol
      albsat
      albdry

      pct_sand
      pct_clay
      organic
      organic_max - scalar

    snicar_data radiation params:
      ss_alb_oc1, 
      asm_prm_oc1,
      ext_cff_mss_oc1,
      ss_alb_oc2,
      asm_prm_oc2,
      ext_cff_mss_oc2,
      ss_alb_dst1, 
      asm_prm_dst1, 
      ext_cff_mss_dst1
      ss_alb_dst2, 
      asm_prm_dst2, 
      ext_cff_mss_dst2
      ss_alb_dst3, 
      asm_prm_dst3,
      ext_cff_mss_dst3
      ss_alb_dst4,
      asm_prm_dst4,
      ext_cff_mss_dst4

      ss_alb_ice_drc
      asm_prm_ice_drc
      ext_cff_mss_ice_drc
      ss_alb_ice_dfs
      asm_prm_ice_dfs
      ext_cff_mss_ice_dfs
      ss_alb_bc_mam
      asm_prm_bc_mam
      ext_cff_mss_bc_mam
      ss_alb_bc_mam
      asm_prm_bc_mam
      ext_cff_mss_bc_mam

      bcenh

      snow radius lookup table
      snowage_tau
      snowage_kappa
      snowage_drdt0


    PFT data:
      fnr;   
      act25
      kcha;  
      koha;  
      cpha;  
      vcmaxha
      jmaxha;
      tpuha; 
      lmrha; 
      vcmaxhd
      jmaxhd;   
      tpuhd;    
      lmrhd;    
      lmrse;    
      qe;       
      theta_cj; 
      bbbopt;   
      mbbopt;   
      c3psn;    
      slatop;   
      leafcn;   
      flnr;     
      fnitr;    
      dleaf;    
      smpso;    
      smpsc;    
      tc_stress;
      z0mr;     
      displar;  
      xl;       
      roota_par;
      rootb_par;
      rholvis;  
      rholnir;  
      rhosvis;  
      rhosnir;  
      taulvis;  
      taulnir;  
      tausvis;  
      tausnir;  


  */
}

