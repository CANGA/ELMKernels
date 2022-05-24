
/*
Temperature is represented as a pentadiagonal system

5-band
matrix (diag-2, diag-1, diag, diag+1, diag+2)

Non-zero pattern of bmatrix:
       SNOW-LAYERS
           |
           |  STANDING-SURFACE-WATER
           |         |
           |         |              SOIL-LAYERS
           |         |                  |
           v         v                  v
     -5 -4 -3 -2 -1| 0| 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
     ==============================================================
 -5 | x  x         |  |                                            |
 -4 | x  x  x      |  |                                            |
 -3 |    x  x  x   |  |                                            |
 -2 |       x  x  x|  |                                            |
 -1 |          x  x|  | x                                          |
     ==============================================================
  0 |              | x| x                                          |
     ==============================================================
  1 |             x| x| x  x                                       |
  2 |              |  | x  x  x                                    |
  3 |              |  |    x  x  x                                 |
  4 |              |  |       x  x  x                              |
  5 |              |  |          x  x  x                           |
  6 |              |  |             x  x  x                        |
  7 |              |  |                x  x  x                     |
  8 |              |  |                   x  x  x                  |
  9 |              |  |                      x  x  x               |
 10 |              |  |                         x  x  x            |
 11 |              |  |                            x  x  x         |
 12 |              |  |                               x  x  x      |
 13 |              |  |                                  x  x  x   |
 14 |              |  |                                     x  x  x|
 15 |              |  |                                        x  x|
     ==============================================================


INDEX CONVERSION::

ELM      this model

         bmatrix    bmatrix_snow   bmatrix_soil   bmatrix_ssw
 i
-5 snow  0          0              --             -- 
-4  |    1          1              --             -- 
-3  |    2          2              --             -- 
-2  |    3          3              --             -- 
-1  |    4          4              --             -- 
 0 ssw   5          --             --             0  
 1  |    6          --             0              -- 
 2  |    7          --             1              -- 
 3  |    8          --             2              -- 
 4  |    9          --             3              -- 
 5  |    10         --             4              -- 
 6  |    11         --             5              -- 
 7  |    12         --             6              -- 
 8  |    13         --             7              -- 
 9  |    14         --             8              -- 
10  |    15         --             9              -- 
11  |    16         --             10             -- 
12  |    17         --             11             -- 
13  |    18         --             12             -- 
14  |    19         --             13             -- 
15 soil  20         --             14             -- 

conversion from ELM
 i       (i+5)      (i+5)          (i-1)          --
conversion from bmatrix
         b          (b)          (b-6)            --

main matrix:
in ELM: bmatrix(ncells, nband, nlevgrnd+nlevsno+1) banded matrix for numerical solution of temperature
here: bmatrix(ncells, nlevgrnd+nlevsno+1, nband) banded matrix for numerical solution of temperature

submatrix original ELM bounds:
bmatrix_snow(bounds%begc:bounds%endc,nband,-nlevsno:-1      )  ! block-diagonal matrix for snow layers
bmatrix_ssw(bounds%begc:bounds%endc,nband,       0:0       )   ! block-diagonal matrix for standing surface water
bmatrix_soil(bounds%begc:bounds%endc,nband,       1:nlevgrnd)  ! block-diagonal matrix for soil layers
bmatrix_snow_soil(bounds%begc:bounds%endc,nband,-1:-1)         ! off-diagonal matrix for snow-soil interaction
bmatrix_ssw_soil(bounds%begc:bounds%endc,nband, 0:0 )          ! off-diagonal matrix for standing surface water-soil interaction
bmatrix_soil_snow(bounds%begc:bounds%endc,nband, 1:1 )         ! off-diagonal matrix for soil-snow interaction
bmatrix_soil_ssw(bounds%begc:bounds%endc,nband, 1:1 )          ! off-diagonal matrix for soil-standing surface water interaction
*/

#pragma once

#include "helper_functions.hh"

namespace ELM::soil_temp_lhs {

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3>
  void set_LHS(const double& dtime,
               const ArrayI1 snl,
               const ArrayD1 dz_h2osfc,
               const ArrayD1 c_h2osfc,
               const ArrayD1 tk_h2osfc,
               const ArrayD1 frac_h2osfc,
               const ArrayD1 frac_sno_eff,
               const ArrayD1 dhsdT,
               const ArrayD2 z,
               const ArrayD2 fact,
               const ArrayD2 tk,
               ArrayD3 lhs_matrix)
  {
    using ELMdims::nlevsno;
    using ELMdims::nlevgrnd;
    using ELMdims::nband;
    using ELM::Utils::create;

    auto bmatrix_snow = create<ArrayD3>("bmatrix_snow", snl.extent(0), nlevsno, nband);
    auto bmatrix_soil = create<ArrayD3>("bmatrix_soil", snl.extent(0), nlevgrnd, nband);
    auto bmatrix_ssw = create<ArrayD2>("bmatrix_ssw", snl.extent(0), nband);
    auto bmatrix_snow_soil = create<ArrayD2>("bmatrix_snow_soil", snl.extent(0), nband);
    auto bmatrix_ssw_soil = create<ArrayD2>("bmatrix_ssw_soil", snl.extent(0), nband);
    auto bmatrix_soil_snow = create<ArrayD2>("bmatrix_soil_snow", snl.extent(0), nband);
    auto bmatrix_soil_ssw = create<ArrayD2>("bmatrix_soil_ssw", snl.extent(0), nband);

    auto kernel = [=] (const int& c) {

      detail::get_matrix_snow(c, snl, dhsdT, z, fact, tk, bmatrix_snow);
      detail::get_matrix_snow_soil(c, snl, z, fact, tk, bmatrix_snow_soil);
      detail::get_matrix_soil(c, snl, dhsdT, frac_sno_eff, frac_h2osfc, dz_h2osfc,
          tk_h2osfc, z, fact, tk, bmatrix_soil);
      detail::get_matrix_soil_snow(c, snl, frac_sno_eff, z, fact, tk, bmatrix_soil_snow);
      detail::get_matrix_ssw(c, dtime, dz_h2osfc, c_h2osfc, tk_h2osfc, dhsdT, z,
          bmatrix_ssw);
      detail::get_matrix_ssw_soil(c, dtime, dz_h2osfc, c_h2osfc, tk_h2osfc, z,
          bmatrix_ssw_soil);
      detail::get_matrix_soil_ssw(c, dtime, frac_h2osfc, dz_h2osfc, tk_h2osfc, fact, z,
          bmatrix_soil_ssw);
      detail::assemble_lhs(c, bmatrix_snow_soil, bmatrix_ssw_soil, bmatrix_soil_snow,
          bmatrix_soil_ssw, bmatrix_ssw, bmatrix_snow, bmatrix_soil, lhs_matrix);
    };

    invoke_kernel(kernel, std::make_tuple(snl.extent(0)), "soil_temp_lhs::set_LHS");
  }

} // namespace ELM::soil_temp_lhs


namespace ELM::soil_temp_lhs::detail {

  static constexpr double cnfac{0.5}; // Crank Nicholson factor between 0 and 1

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void get_matrix_snow(const int& c,
                       const ArrayI1 snl,
                       const ArrayD1 dhsdT,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD3 bmatrix_snow)
  {
    using ELMdims::nlevsno;

    for (int i = 0; i < bmatrix_snow.extent(1); ++i) {
      for (int j = 0; j < bmatrix_snow.extent(2); ++j) {
        bmatrix_snow(c, i, j) = 0.0;
      } 
    }

    if (snl(c) > 0) {
      const int top = nlevsno - snl(c);
      double dzp = z(c, top + 1) - z(c, top);
      bmatrix_snow(c, top, 3) = 0.0;
      bmatrix_snow(c, top, 2) = 1.0 + (1.0 - cnfac) * fact(c, top) * tk(c, top) / dzp - fact(c, top) * dhsdT(c);
      if (snl(c) > 1) {
        bmatrix_snow(c, top, 1) = - (1.0 - cnfac) * fact(c, top) * tk(c, top) / dzp;
      }

      for (int i = top + 1; i < nlevsno; ++i) {
        double dzm = z(c, i) - z(c, i - 1);
        dzp = z(c, i + 1) - z(c, i);
        bmatrix_snow(c, i, 3) = - (1.0 - cnfac) * fact(c, i) * tk(c, i - 1) / dzm;
        bmatrix_snow(c, i, 2) = 1.0 + (1.0 - cnfac) * fact(c, i) * (tk(c, i) / dzp + tk(c, i - 1) / dzm);
        if (i != nlevsno - 1) {
          bmatrix_snow(c, i, 1) = - (1.0 - cnfac) * fact(c, i) * tk(c, i) / dzp;
        }
      }
    }
  }


  template <typename ArrayI1, typename ArrayD2>
  ACCELERATE
  void get_matrix_snow_soil(const int& c,
                       const ArrayI1 snl,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD2 bmatrix_snow_soil)
  {
    using ELMdims::nlevsno;

    for (int i = 0; i < bmatrix_snow_soil.extent(1); ++i) {
      bmatrix_snow_soil(c, i) = 0.0;
    }

    if (snl(c) > 0) {
      bmatrix_snow_soil(c, 0) = - (1.0 - cnfac) * fact(c, nlevsno - 1) * tk(c, nlevsno - 1) /
        (z(c, nlevsno) - z(c, nlevsno-1));
    }
  }


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void get_matrix_soil(const int& c,
                       const ArrayI1 snl,
                       const ArrayD1 dhsdT,
                       const ArrayD1 frac_sno_eff,
                       const ArrayD1 frac_h2osfc,
                       const ArrayD1 dz_h2osfc,
                       const ArrayD1 tk_h2osfc,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD3 bmatrix_soil)
  {
    using ELMdims::nlevsno;
    using ELMdims::nlevgrnd;

    for (int i = 0; i < bmatrix_soil.extent(1); ++i) {
      for (int j = 0; j < bmatrix_soil.extent(2); ++j) {
        bmatrix_soil(c, i, j) = 0.0;
      } 
    }

    if (snl(c) == 0) {
      double dzp = z(c, nlevsno + 1) - z(c, nlevsno);
      bmatrix_soil(c, 0, 2) = 1.0 + (1.0 - cnfac) * fact(c, nlevsno) * tk(c, nlevsno) / dzp - fact(c, nlevsno) * dhsdT(c);
      bmatrix_soil(c, 0, 1) = - (1.0 - cnfac) * fact(c, nlevsno) * tk(c, nlevsno) / dzp;
    } else {
      // this is the snow/soil interface layer
      double dzm = z(c, nlevsno) - z(c, nlevsno - 1);
      double dzp = z(c, nlevsno + 1) - z(c, nlevsno);
      bmatrix_soil(c, 0, 2) = 1.0 + (1.0 - cnfac) * fact(c,nlevsno) * (tk(c, nlevsno) /
        dzp + frac_sno_eff(c) * tk(c, nlevsno - 1) / dzm) - (1.0 - frac_sno_eff(c)) *
        fact(c, nlevsno) * dhsdT(c);
      bmatrix_soil(c, 0, 1) = - (1.0 - cnfac) * fact(c, nlevsno) * tk(c, nlevsno) / dzp;
    }

    for (int i = 1; i < nlevgrnd - 1; ++i) {
      int offset = i + nlevsno;
      double dzm = z(c, offset) - z(c, offset - 1);
      double dzp = z(c, offset + 1) - z(c, offset);
      bmatrix_soil(c, i, 3) = - (1.0 - cnfac) * fact(c, offset) * tk(c, offset - 1) / dzm;
      bmatrix_soil(c, i, 2) = 1.0 + (1.0 - cnfac) * fact(c, offset) *
        (tk(c, offset) / dzp + tk(c, offset - 1) / dzm);
      bmatrix_soil(c, i, 1) = - (1.0 - cnfac) * fact(c, offset) * tk(c, offset) / dzp;
    }

    const int bot = nlevgrnd + nlevsno - 1;
    double dzm = z(c, bot) - z(c, bot - 1);
    bmatrix_soil(c, nlevgrnd - 1, 3) = - (1.0 - cnfac) * fact(c, bot) * tk(c, bot - 1) / dzm;
    bmatrix_soil(c, nlevgrnd - 1, 2) = 1.0 + (1.0 - cnfac) * fact(c, bot) * tk(c, bot - 1) / dzm;
    bmatrix_soil(c, nlevgrnd - 1, 1) = 0.0;

    // diagonal element correction for presence of h2osfc
    if (frac_h2osfc(c) != 0.0) {
      dzm = 0.5 * dz_h2osfc(c) + z(c, nlevsno);
      bmatrix_soil(c, 0, 2) += frac_h2osfc(c) *
        ((1.0 - cnfac) * fact(c, nlevsno) * tk_h2osfc(c) / dzm + fact(c, nlevsno) * dhsdT(c));
    }
  }


  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_soil_snow(const int& c,
                       const ArrayI1 snl,
                       const ArrayD1 frac_sno_eff,
                       const ArrayD2 z,
                       const ArrayD2 fact,
                       const ArrayD2 tk,
                       ArrayD2 bmatrix_soil_snow)
  {
    using ELMdims::nlevsno;

    for (int i = 0; i < bmatrix_soil_snow.extent(1); ++i) {
      bmatrix_soil_snow(c, i) = 0.0;
    }

    if (snl(c) == 0) {
      bmatrix_soil_snow(c, 4) = 0.0;
    } else {
      // this is the snow/soil interface layer
      double dzm = (z(c, nlevsno) - z(c, nlevsno - 1));
      bmatrix_soil_snow(c, 4) = -frac_sno_eff(c) * (1.0 - cnfac) * fact(c, nlevsno) * tk(c,nlevsno - 1) / dzm;
    }
  }


  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_ssw(const int& c,
                      const double& dtime,
                      const ArrayD1 dz_h2osfc,
                      const ArrayD1 c_h2osfc,
                      const ArrayD1 tk_h2osfc,
                      const ArrayD1 dhsdT,
                      const ArrayD2 z,
                      ArrayD2 bmatrix_ssw)
  {
    for (int i = 0; i < bmatrix_ssw.extent(1); ++i) {
      bmatrix_ssw(c, i) = 0.0;
    }

    bmatrix_ssw(c, 2) = 1.0 + (1.0 - cnfac) * (dtime / c_h2osfc(c)) * tk_h2osfc(c) /
        (0.5 * dz_h2osfc(c) + z(c, nlevsno)) - (dtime / c_h2osfc(c)) * dhsdT(c);
  }


  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_ssw_soil(const int& c,
                           const double& dtime,
                           const ArrayD1 dz_h2osfc,
                           const ArrayD1 c_h2osfc,
                           const ArrayD1 tk_h2osfc,
                           const ArrayD2 z,
                           ArrayD2 bmatrix_ssw_soil)
  {
    for (int i = 0; i < bmatrix_ssw_soil.extent(1); ++i) {
      bmatrix_ssw_soil(c, i) = 0.0;
    }

    bmatrix_ssw_soil(c, 1) = - (1.0 - cnfac) * (dtime / c_h2osfc(c)) * tk_h2osfc(c) /
        (0.5 * dz_h2osfc(c) + z(c, nlevsno));
  }


  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void get_matrix_soil_ssw(const int& c,
                           const double& dtime,
                           const ArrayD1 frac_h2osfc,
                           const ArrayD1 dz_h2osfc,
                           const ArrayD1 tk_h2osfc,
                           const ArrayD2 fact,
                           const ArrayD2 z,
                           ArrayD2 bmatrix_soil_ssw)
  {
    for (int i = 0; i < bmatrix_soil_ssw.extent(1); ++i) {
      bmatrix_soil_ssw(c, i) = 0.0;
    }

    if (frac_h2osfc(c) != 0.0) {
      bmatrix_soil_ssw(c, 3) = -frac_h2osfc(c) * (1.0 - cnfac) * fact(c, nlevsno) *
        tk_h2osfc(c) / (0.5 * dz_h2osfc(c) + z(c, nlevsno));
    }
  }


  template <typename ArrayD2, typename ArrayD3>
  ACCELERATE
  void assemble_lhs(const int& c,
                    const ArrayD2 bmatrix_snow_soil,
                    const ArrayD2 bmatrix_ssw_soil,
                    const ArrayD2 bmatrix_soil_snow,
                    const ArrayD2 bmatrix_soil_ssw,
                    const ArrayD2 bmatrix_ssw,
                    const ArrayD3 bmatrix_snow,
                    const ArrayD3 bmatrix_soil,
                    ArrayD3 lhs_matrix)
  {
    using ELMdims::nlevsno;
    using ELMdims::nlevgrnd;

    // zero out matrix
    for (int i = 0; i < lhs_matrix.extent(1); ++i) {
      for (int j = 0; j < lhs_matrix.extent(2); ++j) {
        lhs_matrix(c, i, j) = 0.0;
      } 
    }

    // SNOW
    // bmatrix(c,2:3,-5   ) = bmatrix_snow(c,2:3,-5   )
    // for lev idx = 0 and band idx 1 and 2
    for (int bnd : {1, 2})
      lhs_matrix(c, 0, bnd) = bmatrix_snow(c, 0, bnd);

    // bmatrix(c,2:4,-4:-2) = bmatrix_snow(c,2:4,-4:-2)
    // for lev idx 1 - 3 and band idx 1 - 3
    for (int lev : {1, 2, 3}) {
      for (int bnd : {1, 2, 3}) {
        lhs_matrix(c, lev, bnd) = bmatrix_snow(c, lev, bnd);
      }
    }

    // bmatrix(c,3:4,-1   ) = bmatrix_snow(c,3:4,-1   )
    // for lev idx 4 and band idx 2 - 3
    for (int bnd : {2, 3}) {
      lhs_matrix(c, 4, bnd) = bmatrix_snow(c, 4, bnd);
    }

    // SNOW-SOIL
    // bmatrix(c,1,-1) = bmatrix_snow_soil(c,1,-1)
    // for lev idx 4 and bnd idx 0
    lhs_matrix(c, 4, 0) = bmatrix_snow_soil(c, 0);


    // SSW
    // bmatrix(c,3,0) = bmatrix_ssw(c,3,0)
    // for lev idx 5 and bnd idx 2
    lhs_matrix(c, 5, 2) = bmatrix_ssw(c, 2);

    // SSW-SOIL
    // bmatrix(c,2,0) = bmatrix_ssw_soil(c,2,0)
    // for lev idx 5 and bnd idx 1
    lhs_matrix(c, 5, 1) = bmatrix_ssw_soil(c, 1);

    // SOIL
    // bmatrix(c,2:3,1           )  = bmatrix_soil(c,2:3,1           )
    // for 
    // lhs_matrix lev idx 6
    // bmatrix_soil lev idx 0
    // and band idx 1 - 2
    for (int bnd : {1, 2}) {
      lhs_matrix(c, 6, bnd) = bmatrix_soil(c, 0, bnd);
    }

    // bmatrix(c,2:4,2:nlevgrnd-1)  = bmatrix_soil(c,2:4,2:nlevgrnd-1)
    // for 
    // lhs_matrix lev idx 7 - 19
    // bmatrix_snow lev idx 1 - 13
    // and band idx 1 - 3
    for (int lev = nlevsno + 2; lev < nlevgrnd + nlevsno; ++lev) {
      for (int bnd : {1, 2, 3}) {
        lhs_matrix(c, lev, bnd) = bmatrix_soil(c, lev - 6, bnd);
      }
    }

    // bmatrix(c,3:4,nlevgrnd    )  = bmatrix_soil(c,3:4,nlevgrnd    )
    // for 
    // lhs_matrix lev idx 20
    // bmatrix_soil lev idx 14
    // and band idx 2 - 3
    for (int bnd : {2, 3}) {
      lhs_matrix(c, nlevsno + nlevgrnd, bnd) = bmatrix_soil(c, nlevsno + nlevgrnd - 6, bnd);
    }

    // SOIL-SNOW
    // bmatrix(c,5,1)  = bmatrix_soil_snow(c,5,1)
    // for lev idx 6 and bnd idx 4
    lhs_matrix(c, nlevsno + 1, 4) = bmatrix_soil_snow(c, 4);


    // SOIL-SSW
    // bmatrix(c,4,1)  = bmatrix_soil_ssw(c,4,1)
    // for lev idx 6 and bnd idx 3
    lhs_matrix(c, nlevsno + 1, 3) = bmatrix_soil_ssw(c, 3);
  }






  // basic implementation of a pentadiagonal linear solver
  // will perform well on CPU, possibly poorly on GPU
  // should implement cyclic reduction to maximize GPU performance
  // S. S. Askar, A. A. Karawia, "On Solving Pentadiagonal Linear Systems via Transformations",
  // Mathematical Problems in Engineering, vol. 2015, Article ID 232456, 9 pages, 2015.
  // https://doi.org/10.1155/2015/232456
  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void solve(const int& snl,
             const ArrayD2 lhs_matrix,
             const ArrayD1 rhs_vector)
  {
    using ELMdims::nlevsno;
    using ELMdims::nlevgrnd;
    using ELM::Utils::create;

    // maximum size of system
    const int N = nlevgrnd + nlevsno + 1;
    // top active layer
    const int top = nlevsno - snl;

    // working vectors
    // temporarily allocated here for simplicity
    auto A = create<ArrayD1>("A", N - 1); // A(Ai, ..., An-1)
    auto B = create<ArrayD1>("B", N - 2); // B(Bi, ..., Bn-2)
    auto Z = create<ArrayD1>("Z", N); // Z(Zi, ..., Zn)
    auto Y = create<ArrayD1>("Y", N - 1); // Y(Yi+1, ..., Yn)
    auto U = create<ArrayD1>("U", N); // U(Ui, ..., Un)
    // solution vector
    auto X = create<ArrayD1>("X", N);

    // form coefficients
    // forward sweep
    // top layer calculations
    U(top) = lhs_matrix(top, 2);
    A(top) = lhs_matrix(top, 1) / U(top);
    B(top) = lhs_matrix(top, 0) / U(top);
    Z(top) = rhs_vector(top) / U(top);

    // 2nd from top
    Y(top) = lhs_matrix(top + 1, 3);
    U(top + 1) = lhs_matrix(top + 1, 2) - A(top) * Y(top);
    A(top + 1) = (lhs_matrix(top + 1, 1) - B(top) * Y(top)) / U(top + 1);
    B(top + 1) = lhs_matrix(top + 1, 0) / U(top + 1);
    Z(top + 1) = (rhs_vector(top + 1) - Z(top) * Y(top)) / U(top + 1);

    // for cells (top + 2 : N - 3)
    for (int i = top + 2; i < N - 2; ++i) {
      Y(i - 1) = lhs_matrix(i, 3) - A(i - 2) * lhs_matrix(i, 4);
      U(i) = lhs_matrix(i, 2) - B(i - 2) * lhs_matrix(i, 4) - A(i - 1) * Y(i - 1);
      A(i) = (lhs_matrix(i, 1) - B(i - 1) * Y(i - 1)) / U(i);
      B(i) = lhs_matrix(i, 0) / U(i);
      Z(i) = (rhs_vector(i) - Z(i - 2) * lhs_matrix(i, 4) - Z(i - 1) * Y(i - 1)) / U(i);
    }

    // 2nd from bottom
    Y(N - 3) = lhs_matrix(N - 2, 3) - A(N - 4) * lhs_matrix(N - 2, 4);
    U(N - 2) = lhs_matrix(N - 2, 2) - B(N - 4) * lhs_matrix(N - 2, 4) -
      A(N - 3) * Y(N - 3);
    A(N - 2) = (lhs_matrix(N - 2, 1) - B(N - 3) * Y(N - 3)) / U(N - 2);

    // bottom cell
    Y(N - 2) = lhs_matrix(N - 1, 3) - A(N - 3) * lhs_matrix(N - 1, 4);
    U(N - 1) = lhs_matrix(N - 1, 2) - B(N - 3) * lhs_matrix(N - 1, 4) -
      A(N - 2) * Y(N - 2);
    Z(N - 2) = (rhs_vector(N - 2) - Z(N - 3) * lhs_matrix(N - 2, 4) -
      Z(N - 3) * Y(N - 3)) / U(N - 2);
    Z(N - 1) = (rhs_vector(N - 1) - Z(N - 2) * lhs_matrix(N - 1, 4) -
      Z(N - 2) * Y(N - 2)) / U(N - 1);

    // obtain solution via backward substitution
    X(N - 1) = Z(N - 1);
    X(N - 2) = Z(N - 2) - A(N - 2) * X(N - 1);
    for (int i = N - 3; i >= 0; --i)
      X(i) = Z(i) - A(i) * X(i + 1) - B(i) * X(i + 2);

  }


} // namespace ELM::soil_temp_lhs::detail
