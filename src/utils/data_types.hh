
#pragma once

#include "pft_data.h"
#include "elm_state.h"

//#ifndef ENABLE_KOKKOS
//
//namespace ELM { template <typename T, size_t D> class Array; }
//typedef ELM::Array<bool, 1> ViewB1;
//typedef ELM::Array<int, 1> ViewI1;
//typedef ELM::Array<int, 2> ViewI2;
//typedef ELM::Array<std::string, 1> ViewS1;
//typedef ELM::Array<double, 1> ViewD1;
//typedef ELM::Array<double, 2> ViewD2;
//typedef ELM::Array<double, 3> ViewD3;
//typedef ViewI1 h_ViewI1;
//typedef ViewD1 h_ViewD1;
//typedef ViewD2 h_ViewD2;
//typedef ViewD3 h_ViewD3;
//typedef ELM::Array<ELM::PFTDataPSN, 1> ViewPSN1;
//
//typedef ELM::ELMState<ViewB1, ViewI1, ViewI2, ViewD1, ViewD2, ViewD3, ViewPSN1> ELMStateType;
//
//#else

typedef Kokkos::View<bool *> ViewB1;
typedef Kokkos::View<int *> ViewI1;
typedef Kokkos::View<int **> ViewI2;
typedef Kokkos::View<std::string *> ViewS1;
typedef Kokkos::View<double *> ViewD1;
typedef Kokkos::View<double **> ViewD2;
typedef Kokkos::View<double ***> ViewD3;
typedef ViewI1::HostMirror h_ViewI1;
typedef ViewD1::HostMirror h_ViewD1;
typedef ViewD2::HostMirror h_ViewD2;
typedef ViewD3::HostMirror h_ViewD3;
typedef Kokkos::View<ELM::PFTDataPSN *> ViewPSN1;

typedef ELM::ELMState<ViewB1, ViewI1, ViewI2, ViewD1, ViewD2, ViewD3, ViewPSN1> ELMStateType;

//#endif
