#pragma once



namespace ELM::LND {
// from landunit_varcon.F90
static constexpr int istsoil{1};
static constexpr int istcrop{2};
static constexpr int istice{3};
static constexpr int istice_mec{4};
static constexpr int istdlak{5};
static constexpr int istwet{6};
static constexpr int isturb_MIN{7};
static constexpr int isturb_tbd{7};
static constexpr int isturb_hd{8};
static constexpr int isturb_md{9};
static constexpr int isturb_MAX{9};
static constexpr int max_lunit{9};
static constexpr int landunit_name_length{40};
static constexpr int numurbl{isturb_MAX - isturb_MIN + 1};

// from column_varcon.F90
static constexpr int icol_roof{isturb_MIN * 10 + 1};
static constexpr int icol_sunwall{isturb_MIN * 10 + 2};
static constexpr int icol_shadewall{isturb_MIN * 10 + 3};
static constexpr int icol_road_imperv{isturb_MIN * 10 + 4};
static constexpr int icol_road_perv{isturb_MIN * 10 + 5};

} // namespace ELM::LND


// simple struct to hold landtype data
// currently used as a placeholder
namespace ELM {

struct LandType {

  LandType() : ltype{1}, ctype{0}, vtype{2}, urbpoi{false}, lakpoi{false} {}

  // land unit, urban unit, vegetation unit
  int ltype, ctype, vtype;
  // urban point, lake point
  bool urbpoi, lakpoi;
};

} // namespace ELM
