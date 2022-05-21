

#pragma once

namespace ELM::soil_matrix {

template <class Array_t> Array_t create(const std::string &name, int D0)
{ return Array_t(name, D0); }
template <class Array_t> Array_t create(const std::string &name, int D0, int D1)
{ return Array_t(name, D0, D1); }

/* DESCRIPTION:
Setup the RHS-Vector for the numerical solution of temperature for snow,
standing surface water and soil layers.
          |===========|
          |   Snow    |
          !===========|
rvector = |    SSW    |
          !===========|
          !   Soil    |
          !===========|

 rt(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)        ! "r" vector for tridiagonal solution
 fn_h2osfc(bounds%begc:bounds%endc)                      ! heat diffusion through standing-water/soil interface [W/m2]
 rt_snow(bounds%begc:bounds%endc,-nlevsno:-1)            ! RHS vector corresponding to snow layers
 rt_ssw(bounds%begc:bounds%endc,1)                       ! RHS vector corresponding to standing surface water
 rt_soil(bounds%begc:bounds%endc,1:nlevgrnd)             ! RHS vector corresponding to soil layer

 // for h2osfc
 const ArrayD1 tk_h2osfc,
 const ArrayD1 t_h2osfc,
 const ArrayD1 dz_h2osfc,
 const ArrayD1 c_h2osfc,
 const ArrayD1 hs_h2osfc,

 // later
 ArrayD2 rvector
*/
  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  void set_rhs(const double& dtime,
               const ArrayI1 snl,
               const ArrayD1 hs_top_snow,
               const ArrayD1 dhsdT,
               const ArrayD1 hs_soil,
               const ArrayD1 frac_sno_eff,
               const ArrayD2 t_soisno,
               const ArrayD2 fact,
               const ArrayD2 fn,
               const ArrayD2 sabg_lyr,
               const ArrayD2 z)
  {
    using ELMdims::nlevsno;
    

    auto fn_h2osfc = create<ArrayD1>("fn_h2osfc", snl.extent(0));
    auto rt_snow = create<ArrayD2>("rt_snow", snl.extent(0), nlevsno);
    auto rt_ssw = create<ArrayD1>("rt_ssw", snl.extent(0));
    auto rt_soil = create<ArrayD2>("rt_soil", snl.extent(0), nlevgrnd);

    auto kernel = [=] (const int& c) {
      set_rhs_snow(c, snl, hs_top_snow, dhsdT, t_soisno,
                   fact, fn, sabg_lyr, rt_snow);
      //set_rhs_ssw();
      set_rhs_soil(c, snl, hs_soil, hs_top_snow, frac_sno_eff,
                   dhsdT, t_soisno, fact, fn, sabg_lyr,
                   rt_soil);
    };

    invoke_kernel(kernel, std::make_tuple(snl.extent(0)), "set_rhs");

    for (int i = 0; i < ELMdims::nlevgrnd; ++i)
      std::cout << "set_rhs: " << i << "  " << rt_soil(0,i) << std::endl; 
    
  }

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void set_rhs_snow(const int& c,
                    const ArrayI1 snl,
                    const ArrayD1 hs_top_snow,
                    const ArrayD1 dhsdT,
                    const ArrayD2 t_soisno,
                    const ArrayD2 fact,
                    const ArrayD2 fn,
                    const ArrayD2 sabg_lyr,
                    ArrayD2 rt_snow)
  {
    using ELMdims::nlevsno;

    const int top = nlevsno - snl(c);
    if (top > 0) {
      
      // set inactive layers to zero
      for (int i = 0; i < top; ++i) {
        rt_snow(c, i) = 0.0;
      }

      // set top active layer
      rt_snow(c, top) = t_soisno(c, top) + fact(c, top) * (hs_top_snow(c) -
          dhsdT(c) * t_soisno(c, top) + cnfac * fn(c, top));

      for (int i = top + 1; i < nlevsno; ++i) {
        rt_snow(c, i) = t_soisno(c, i) + cnfac * fact(c, i) * (fn(c, i) - fn(c, i-1)) +
            fact(c, i) * sabg_lyr(c, i);
      }
    } else {
      // set inactive layers to zero
      for (int i = 0; i < nlevsno; ++i) {
        rt_snow(c, i) = 0.0;
      }
    }
  }

  // rhs for h2osfc
  // surface water layer has two coefficients
  template <typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void set_rhs_ssw(const int& c,
                   const double& dtime,
                   const ArrayD1 tk_h2osfc,
                   const ArrayD1 t_h2osfc,
                   const ArrayD1 dz_h2osfc,
                   const ArrayD1 c_h2osfc,
                   const ArrayD1 hs_h2osfc,
                   const ArrayD1 dhsdT,
                   const ArrayD2 t_soisno,
                   const ArrayD2 z,
                   ArrayD1 fn_h2osfc,
                   ArrayD1 rt_ssw)
  {
    using ELMdims::nlevsno;
    
    fn_h2osfc(c) = tk_h2osfc(c) * (t_soisno(c, nlevsno) - t_h2osfc(c)) /
        (0.5 * dz_h2osfc(c) + z(c, nlevsno));
    rt_ssw(c) = t_h2osfc(c) + (dtime / c_h2osfc(c)) * (hs_h2osfc(c) - dhsdT(c) *
        t_h2osfc(c) + cnfac * fn_h2osfc(c));

  }

  template <typename ArrayI1, typename ArrayD1, typename ArrayD2>
  ACCELERATE
  void set_rhs_soil(const int& c,
                    const ArrayI1 snl,
                    const ArrayD1 hs_soil,
                    const ArrayD1 hs_top_snow,
                    const ArrayD1 frac_sno_eff,
                    const ArrayD1 dhsdT,
                    const ArrayD2 t_soisno,
                    const ArrayD2 fact,
                    const ArrayD2 fn,
                    const ArrayD2 sabg_lyr,
                    ArrayD2 rt_soil)
  {
    using ELMdims::nlevsno;
    using ELMdims::nlevgrnd;

    // top subsurface layer
    if (snl(c) == 0) {
      rt_soil(c, 0) = t_soisno(c, nlevsno) + fact(c, nlevsno) * (hs_top_snow(c) -
          dhsdT(c) * t_soisno(c, nlevsno) + cnfac * fn(c,nlevsno));
    } else {
      // this is the snow/soil interface layer
      rt_soil(c, 0) = t_soisno(c, nlevsno) + fact(c, nlevsno) * ((1.0 - frac_sno_eff(c)) *
          (hs_soil(c) - dhsdT(c) * t_soisno(c, nlevsno)) + cnfac * (fn(c, nlevsno) -
          frac_sno_eff(c) * fn(c, nlevsno - 1)));

      rt_soil(c, 0) += frac_sno_eff(c) * fact(c, nlevsno) * sabg_lyr(c, nlevsno);
    }

    // layers from nlevsno + 1 (second from top) to nlevgrnd + nlevsno - 2 (second from bottom)
    const int bot = nlevgrnd + nlevsno - 1;
    for (int j = nlevsno + 1; j < bot; ++j) {
      rt_soil(c,j - nlevsno) = t_soisno(c, j) + cnfac * fact(c, j) * (fn(c, j) - fn(c,j-1));
    }

    // bottom layer
    rt_soil(c, nlevgrnd - 1) = t_soisno(c, bot) - cnfac * fact(c, bot) * fn(c, bot - 1) +
        fact(c, bot) * fn(c, bot);
  }



} // namespace ELM::soil_matrix {
