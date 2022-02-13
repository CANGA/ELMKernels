
#include "snicar_data.h"
#include "snow_snicar.h"


ELM::SnicarData::SnicarData() : 
    ss_alb_oc1("ss_alb_oc1", snow_snicar::numrad_snw),
    asm_prm_oc1("asm_prm_oc1", snow_snicar::numrad_snw),
    ext_cff_mss_oc1("ext_cff_mss_oc1", snow_snicar::numrad_snw),
    ss_alb_oc2("ss_alb_oc2", snow_snicar::numrad_snw),
    asm_prm_oc2("asm_prm_oc2", snow_snicar::numrad_snw),
    ext_cff_mss_oc2("ext_cff_mss_oc2", snow_snicar::numrad_snw),
    ss_alb_dst1("ss_alb_dst1", snow_snicar::numrad_snw),
    asm_prm_dst1("asm_prm_dst1", snow_snicar::numrad_snw),
    ext_cff_mss_dst1("ext_cff_mss_dst1", snow_snicar::numrad_snw),
    ss_alb_dst2("ss_alb_dst2", snow_snicar::numrad_snw),
    asm_prm_dst2("asm_prm_dst2", snow_snicar::numrad_snw),
    ext_cff_mss_dst2("ext_cff_mss_dst2", snow_snicar::numrad_snw),
    ss_alb_dst3("ss_alb_dst3", snow_snicar::numrad_snw),
    asm_prm_dst3("asm_prm_dst3", snow_snicar::numrad_snw),
    ext_cff_mss_dst3("ext_cff_mss_dst3", snow_snicar::numrad_snw),
    ss_alb_dst4("ss_alb_dst4", snow_snicar::numrad_snw),
    asm_prm_dst4("asm_prm_dst4", snow_snicar::numrad_snw),
    ext_cff_mss_dst4("ext_cff_mss_dst4", snow_snicar::numrad_snw),
    ss_alb_snw_drc("ss_alb_snw_drc", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
    asm_prm_snw_drc("asm_prm_snw_drc", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
    ext_cff_mss_snw_drc("ext_cff_mss_snw_drc", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
    ss_alb_snw_dfs("ss_alb_snw_dfs", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
    asm_prm_snw_dfs("asm_prm_snw_dfs", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
    ext_cff_mss_snw_dfs("ext_cff_mss_snw_dfs", snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx),
    ss_alb_bc1("ss_alb_bc1", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw),
    asm_prm_bc1("asm_prm_bc1", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw),
    ext_cff_mss_bc1("ext_cff_mss_bc1", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw),
    ss_alb_bc2("ss_alb_bc2", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw),
    asm_prm_bc2("asm_prm_bc2", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw),
    ext_cff_mss_bc2("ext_cff_mss_bc2", ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw),
    bcenh("bcenh", snow_snicar::idx_bcint_icerds_max+1, snow_snicar::idx_bc_nclrds_max+1, 
      snow_snicar::numrad_snw)  { }


void ELM::read_snicar_data(const Comm_type &comm, const std::string& filename, SnicarData *snicar_data) {

  // get file start idx and size to read for 1D arrays
  {
    std::array<size_t, 1> start = {0};
    std::array<size_t, 1> count = {snow_snicar::numrad_snw};
    Array<double, 1> arr_for_read(snow_snicar::numrad_snw);
  
    // read 1D arrays of size snow_snicar::numrad_snw
    read_and_fill_array(comm, filename, "ss_alb_ocphil", arr_for_read, snicar_data->ss_alb_oc1);
    read_and_fill_array(comm, filename, "asm_prm_ocphil", arr_for_read, snicar_data->asm_prm_oc1);
    read_and_fill_array(comm, filename, "ext_cff_mss_ocphil", arr_for_read, snicar_data->ext_cff_mss_oc1);
    read_and_fill_array(comm, filename, "ss_alb_ocphob", arr_for_read, snicar_data->ss_alb_oc2);
    read_and_fill_array(comm, filename, "asm_prm_ocphob", arr_for_read, snicar_data->asm_prm_oc2);
    read_and_fill_array(comm, filename, "ext_cff_mss_ocphob", arr_for_read, snicar_data->ext_cff_mss_oc2);
    read_and_fill_array(comm, filename, "ss_alb_dust01", arr_for_read, snicar_data->ss_alb_dst1);
    read_and_fill_array(comm, filename, "asm_prm_dust01", arr_for_read, snicar_data->asm_prm_dst1);
    read_and_fill_array(comm, filename, "ext_cff_mss_dust01", arr_for_read, snicar_data->ext_cff_mss_dst1);
    read_and_fill_array(comm, filename, "ss_alb_dust02", arr_for_read, snicar_data->ss_alb_dst2);
    read_and_fill_array(comm, filename, "asm_prm_dust02", arr_for_read, snicar_data->asm_prm_dst2);
    read_and_fill_array(comm, filename, "ext_cff_mss_dust02", arr_for_read, snicar_data->ext_cff_mss_dst2);
    read_and_fill_array(comm, filename, "ss_alb_dust03", arr_for_read, snicar_data->ss_alb_dst3);
    read_and_fill_array(comm, filename, "asm_prm_dust03", arr_for_read, snicar_data->asm_prm_dst3);
    read_and_fill_array(comm, filename, "ext_cff_mss_dust03", arr_for_read, snicar_data->ext_cff_mss_dst3);
    read_and_fill_array(comm, filename, "ss_alb_dust04", arr_for_read, snicar_data->ss_alb_dst4);
    read_and_fill_array(comm, filename, "asm_prm_dust04", arr_for_read, snicar_data->asm_prm_dst4);
    read_and_fill_array(comm, filename, "ext_cff_mss_dust04", arr_for_read, snicar_data->ext_cff_mss_dst4);
  }

  // get file start idx and size to read for 2D arrays
  {
    std::array<size_t, 2> start = {0};
    std::array<size_t, 2> count = {snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx};
    Array<double, 2> arr_for_read(snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx);
  
    // read 2D arrays of size [snow_snicar::numrad_snw, snow_snicar::idx_Mie_snw_mx]
    read_and_fill_array(comm, filename, "ss_alb_ice_drc", arr_for_read, snicar_data->ss_alb_snw_drc);
    read_and_fill_array(comm, filename, "asm_prm_ice_drc", arr_for_read, snicar_data->asm_prm_snw_drc);
    read_and_fill_array(comm, filename, "ext_cff_mss_ice_drc", arr_for_read, snicar_data->ext_cff_mss_snw_drc);
    read_and_fill_array(comm, filename, "ss_alb_ice_dfs", arr_for_read, snicar_data->ss_alb_snw_dfs);
    read_and_fill_array(comm, filename, "asm_prm_ice_dfs", arr_for_read, snicar_data->asm_prm_snw_dfs);
    read_and_fill_array(comm, filename, "ext_cff_mss_ice_dfs", arr_for_read, snicar_data->ext_cff_mss_snw_dfs);
  
  }

  // get file start idx and size to read for 2D arrays
  {
    std::array<size_t, 2> start = {0};
    std::array<size_t, 2> count = {ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw};
    Array<double, 2> arr_for_read(ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw);
  
    // read 2D arrays of size [ELM::snow_snicar::idx_bc_nclrds_max+1, ELM::snow_snicar::numrad_snw]
    read_and_fill_array(comm, filename, "ss_alb_bc_mam", arr_for_read, snicar_data->ss_alb_bc1);
    read_and_fill_array(comm, filename, "asm_prm_bc_mam", arr_for_read, snicar_data->asm_prm_bc1);
    read_and_fill_array(comm, filename, "ext_cff_mss_bc_mam", arr_for_read, snicar_data->ext_cff_mss_bc1);
    read_and_fill_array(comm, filename, "ss_alb_bc_mam", arr_for_read, snicar_data->ss_alb_bc2);
    read_and_fill_array(comm, filename, "asm_prm_bc_mam", arr_for_read, snicar_data->asm_prm_bc2);
    read_and_fill_array(comm, filename, "ext_cff_mss_bc_mam", arr_for_read, snicar_data->ext_cff_mss_bc2);
  }

  // get file start idx and size to read for 3D arrays
  {
    std::array<size_t, 3> start = {0};
    std::array<size_t, 3> count = {snow_snicar::idx_bcint_icerds_max+1, snow_snicar::idx_bc_nclrds_max+1, 
       snow_snicar::numrad_snw};
    Array<double, 3> arr_for_read(snow_snicar::idx_bcint_icerds_max+1, snow_snicar::idx_bc_nclrds_max+1, 
       snow_snicar::numrad_snw);
    // read 3D array
    read_and_fill_array(comm, filename, "bcint_enh_mam", arr_for_read, snicar_data->bcenh);
  }
}

