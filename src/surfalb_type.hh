/*
DESCRIPTION:
Store and initialize variables used in albedo and radiation calculations with reasonable default values.
*/

struct SurfAlb_type {
    
double albgrd_col     = 0.2;
double albgri_col     = 0.2;
double albsod_col     = 0.2;
double albsoi_col     = 0.2;
double albsnd_hst_col = 0.6;
double albsni_hst_col = 0.6;
double albd_patch     = 0.2;
double albi_patch     = 0.2;
double albgrd_pur_col = 0.2;
double albgri_pur_col = 0.2;
double albgrd_bc_col  = 0.2;
double albgri_bc_col  = 0.2;
double albgrd_oc_col  = 0.2;
double albgri_oc_col  = 0.2;
double albgrd_dst_col = 0.2;
double albgri_dst_col = 0.2;
double fabi_patch     = 0.0;
double fabd_patch     = 0.0;
double fabi_sun_patch = 0.0;
double fabd_sun_patch = 0.0;
double fabd_sha_patch = 0.0;
double fabi_sha_patch = 0.0;
double ftdd_patch     = 1.0;
double ftid_patch     = 0.0;
double ftii_patch     = 1.0;
};
