#ifndef PROB_H
#define PROB_H

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

#include "IndexDefines.H"
#include "Constants.H"
#include "EOS.H"
#include "Tagging.H"
#include "Transport.H"
#include "ProblemSpecificFunctions.H"
#include "prob_parm.H"
#include "EB.H"

#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Translation.H>

#include "PeleCUtilities.H"
#include "Geometry.H"

#include <cmath>
#include <string>

void parse_params(ProbParmDevice* prob_parm_device);

class EBAnnularSector
  : public pele::pelec::Geometry::Register<EBAnnularSector>
{
public:
  static std::string identifier()
  {
    return "annular-sector";
  }

  void
  build(
    const amrex::Geometry& geom,
    const int max_coarsening_level) override;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
  const amrex::Real u = 0.0;
  const amrex::Real v = 0.0;
  const amrex::Real w = 0.0;
  const amrex::Real p = prob_parm.p0;
  const amrex::Real T = prob_parm.T0;

  // Species zero initially has mass fraction one.
  amrex::Real massfrac[NUM_SPECIES] = {1.0};

  amrex::Real rho;
  amrex::Real eint;

  auto eos =
    pele::physics::PhysicsType::eos();

  eos.PYT2RE(
    p,
    massfrac,
    T,
    rho,
    eint);

  state(i, j, k, URHO) = rho;
  state(i, j, k, UMX) = rho * u;
  state(i, j, k, UMY) = rho * v;
  state(i, j, k, UMZ) = rho * w;
  state(i, j, k, UEINT) = rho * eint;

  state(i, j, k, UEDEN) =
    rho *
    (
      eint +
      0.5 * (u * u + v * v + w * w)
    );

  state(i, j, k, UTEMP) = T;

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, UFS + n) =
      rho * massfrac[n];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<
    amrex::Real,
    AMREX_SPACEDIM>& turb_fluc)
{
  // ---------------------------------------------------------------------------
  // x-low inlet
  // ---------------------------------------------------------------------------

  if ((idir == 0) && (sgn == 1)) {

    auto eos =
      pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES];

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] =
        s_int[UFS + n] / s_int[URHO];
    }

    amrex::Real p_in = 0.0;

    if (prob_parm.inlet_type == 0) {

      eos.RTY2P(
        s_int[URHO],
        s_int[UTEMP],
        massfrac,
        p_in);

    } else if (prob_parm.inlet_type == 1) {

      p_in = prob_parm.p_inlet;

    } else {

      amrex::Abort(
        "Invalid inlet boundary type specified");
    }

    const amrex::Real T_in =
      prob_parm.T_inlet;

    amrex::Real rho_in;
    amrex::Real eint_in;

    eos.PYT2RE(
      p_in,
      massfrac,
      T_in,
      rho_in,
      eint_in);

    amrex::Real Cs_in;

    eos.RTY2Cs(
      rho_in,
      T_in,
      massfrac,
      Cs_in);

    amrex::Real u_in =
      prob_parm.M_inlet * Cs_in;

    amrex::Real v_in = 0.0;
    amrex::Real w_in = 0.0;

    if (
      prob_parm.turb_inflow_type ==
      TurbInflowType::TurbulentFluctuations) {

      AMREX_D_TERM(
        u_in += turb_fluc[0];,
        v_in += turb_fluc[1];,
        w_in += turb_fluc[2];);

    } else if (
      prob_parm.turb_inflow_type ==
      TurbInflowType::TotalVelocity) {

      AMREX_D_TERM(
        u_in = turb_fluc[0];,
        v_in = turb_fluc[1];,
        w_in = turb_fluc[2];);
    }

    s_ext[URHO] = rho_in;
    s_ext[UMX] = rho_in * u_in;
    s_ext[UMY] = rho_in * v_in;
    s_ext[UMZ] = rho_in * w_in;
    s_ext[UEINT] = rho_in * eint_in;

    s_ext[UEDEN] =
      rho_in *
      (
        eint_in +
        0.5 *
        (
          u_in * u_in +
          v_in * v_in +
          w_in * w_in
        )
      );

    s_ext[UTEMP] = T_in;

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] =
        rho_in * massfrac[n];
    }

  // ---------------------------------------------------------------------------
  // x-high outlet
  // ---------------------------------------------------------------------------

  } else if ((idir == 0) && (sgn == -1)) {

    amrex::Real rho =
      s_int[URHO];

    amrex::Real u =
      s_int[UMX] / rho;

    amrex::Real v =
      s_int[UMY] / rho;

    amrex::Real w =
      s_int[UMZ] / rho;

    amrex::Real eint =
      s_int[UEINT] / rho;

    amrex::Real T =
      s_int[UTEMP];

    amrex::Real massfrac[NUM_SPECIES];

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] =
        s_int[UFS + n] / rho;
    }

    if (prob_parm.outlet_type == 0) {

      // First-order extrapolation.

    } else if (prob_parm.outlet_type == 1) {

      const amrex::Real* prob_hi =
        geomdata.ProbHi();

      const amrex::Real* dx =
        geomdata.CellSize();

      const amrex::Real xd =
        prob_hi[0] - 0.5 * dx[0];

      const amrex::Real xb =
        prob_hi[0];

      const amrex::Real rho_inter =
        rho;

      const amrex::Real u_inter =
        u;

      const amrex::Real v_inter =
        v;

      const amrex::Real w_inter =
        w;

      const amrex::Real T_inter =
        T;

      amrex::Real p_inter = 0.0;
      amrex::Real cs_inter = 0.0;

      auto eos =
        pele::physics::PhysicsType::eos();

      eos.RTY2P(
        rho_inter,
        T_inter,
        massfrac,
        p_inter);

      eos.RTY2Cs(
        rho_inter,
        T_inter,
        massfrac,
        cs_inter);

      const amrex::Real pb =
        prob_parm.p_exit;

      const amrex::Real rhob =
        rho_inter +
        (pb - p_inter) /
        (cs_inter * cs_inter);

      const amrex::Real ub =
        u_inter +
        (p_inter - pb) /
        (rho_inter * cs_inter);

      const amrex::Real vb =
        v_inter;

      const amrex::Real wb =
        w_inter;

      rho =
        (rhob - rho_inter) /
        (xb - xd) *
        (x[0] - xd) +
        rho_inter;

      const amrex::Real p =
        (pb - p_inter) /
        (xb - xd) *
        (x[0] - xd) +
        p_inter;

      eos.RYP2E(
        rho,
        massfrac,
        p,
        eint);

      eos.REY2T(
        rho,
        eint,
        massfrac,
        T);

      u =
        (ub - u_inter) /
        (xb - xd) *
        (x[0] - xd) +
        u_inter;

      v =
        (vb - v_inter) /
        (xb - xd) *
        (x[0] - xd) +
        v_inter;

      w =
        (wb - w_inter) /
        (xb - xd) *
        (x[0] - xd) +
        w_inter;

    } else if (prob_parm.outlet_type == 2) {

      const amrex::Real p =
        prob_parm.p_exit;

      auto eos =
        pele::physics::PhysicsType::eos();

      eos.PYT2RE(
        p,
        massfrac,
        T,
        rho,
        eint);

    } else {

      amrex::Abort(
        "Invalid outlet boundary type specified");
    }

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;

    s_ext[UEDEN] =
      rho *
      (
        eint +
        0.5 * (u * u + v * v + w * w)
      );

    s_ext[UTEMP] = T;

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[UFS + n] =
        rho * massfrac[n];
    }
  }
}

void pc_prob_close();

struct MyProblemSpecificFunctions
  : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  problem_modify_transport_coeffs(
    const bool /*get_xi*/,
    const bool get_mu,
    const bool get_lam,
    const bool /*get_Ddiag*/,
    const bool /*get_chi*/,
    const amrex::Real /*Tloc*/,
    const amrex::Real /*rholoc*/,
    amrex::Real* /*Yloc*/,
    amrex::Real* /*Ddiag*/,
    amrex::Real* /*chi_mix*/,
    amrex::Real& mu,
    amrex::Real& /*xi*/,
    amrex::Real& lam,
    pele::physics::transport::TransParm<
      pele::physics::PhysicsType::eos_type,
      pele::physics::PhysicsType::transport_type> const*
      /*tparm*/,
    ProbParmDevice const& prob_parm,
    const amrex::RealVect& x)
  {
    if (prob_parm.do_sponge_zones) {

      const amrex::Real stretch = 2.0;
      const amrex::Real buffer = 3.0;
      const amrex::Real beg = 0.0;
      const amrex::Real end = prob_parm.L;

      const amrex::Real factor =
        0.5 *
          (
            1.0 +
            std::tanh(
              stretch *
              (beg + buffer - x[0]))
          ) +
        0.5 *
          (
            1.0 +
            std::tanh(
              stretch *
              (x[0] - end + buffer))
          );

      const amrex::Real visc = 1.0;
      const amrex::Real cond = 1.0;

      // Preserves the original problem behavior.
      if (get_mu) {
        mu = factor * visc;
      }

      if (get_lam) {
        lam = factor * cond;
      }
    }
  }
};

using ProblemSpecificFunctions =
  MyProblemSpecificFunctions;

#endif
