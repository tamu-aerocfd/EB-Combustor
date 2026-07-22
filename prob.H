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
#include <cstdio>
#include <string>

void parse_params(ProbParmDevice* prob_parm_device);


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
pc_fuel_species_id(const ProbParmDevice& prob_parm)
{
#if NUM_SPECIES == 1
  return 0;
#else
  return (prob_parm.fuel_species_id >= 0) ? prob_parm.fuel_species_id : H2_ID;
#endif
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
pc_is_fuel_inlet_port(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::GeometryData const& geomdata,
  const ProbParmDevice& prob_parm)
{
  const amrex::Real* dx =
    geomdata.CellSize();

  const amrex::Real dy = x[1] - prob_parm.y_evaporator_center;
  const amrex::Real dz = x[2] - prob_parm.z_evaporator_center;
  const amrex::Real radius_tolerance =
    0.5 * std::sqrt(dx[1] * dx[1] + dx[2] * dx[2]);
  const amrex::Real port_radius =
    prob_parm.r_fuel_line_inner + radius_tolerance;

  return dy * dy + dz * dz <= port_radius * port_radius;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pc_air_massfrac(
  amrex::Real massfrac[NUM_SPECIES],
  const ProbParmDevice& prob_parm)
{
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = 0.0;
  }

#if NUM_SPECIES == 1
  massfrac[0] = 1.0;
#else
  amrex::Real molefrac[NUM_SPECIES] = {0.0};

  const amrex::Real air_mole_frac_sum =
    prob_parm.air_mole_frac_o2 + prob_parm.air_mole_frac_n2;

  molefrac[O2_ID] =
    prob_parm.air_mole_frac_o2 / air_mole_frac_sum;
  molefrac[N2_ID] =
    prob_parm.air_mole_frac_n2 / air_mole_frac_sum;

  auto eos =
    pele::physics::PhysicsType::eos();

  eos.X2Y(
    molefrac,
    massfrac);
#endif
}

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

  amrex::Real massfrac[NUM_SPECIES];
  pc_air_massfrac(
    massfrac,
    prob_parm);

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
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm,
  const amrex::GpuArray<
    amrex::Real,
    AMREX_SPACEDIM>& turb_fluc)
{
  const amrex::Real* prob_lo =
    geomdata.ProbLo();
  const amrex::Real* prob_hi =
    geomdata.ProbHi();
  const amrex::Real x_mid =
    0.5 * (prob_lo[0] + prob_hi[0]);

  const bool is_xlo_boundary =
    (idir == 0) && (x[0] < x_mid);
  const bool is_xhi_boundary =
    (idir == 0) && (x[0] >= x_mid);

  // ---------------------------------------------------------------------------
  // x-low inlet
  // ---------------------------------------------------------------------------

  if (is_xlo_boundary) {

    auto eos =
      pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES];
    pc_air_massfrac(
      massfrac,
      prob_parm);

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

  } else if (is_xhi_boundary) {

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

    auto eos =
      pele::physics::PhysicsType::eos();

    amrex::Real p_inter = 0.0;

    eos.RTY2P(
      rho,
      T,
      massfrac,
      p_inter);

    const bool is_secondary_fuel_inlet =
      pc_is_fuel_inlet_port(x, geomdata, prob_parm);

    if (is_secondary_fuel_inlet) {

      if (prob_parm.debug_fuel_port) {
        printf(
          "Fuel port BC reached at x=(%.16g, %.16g, %.16g), "
          "center=(%.16g, %.16g), inner_radius=%.16g\n",
          static_cast<double>(x[0]),
          static_cast<double>(x[1]),
          static_cast<double>(x[2]),
          static_cast<double>(prob_parm.y_evaporator_center),
          static_cast<double>(prob_parm.z_evaporator_center),
          static_cast<double>(prob_parm.r_fuel_line_inner));
      }

      const int fuel_species_id =
        pc_fuel_species_id(prob_parm);

      const amrex::Real fuel_mass_fraction =
        (prob_parm.mass_frac_fuel < 0.0) ?
          0.0 :
          ((prob_parm.mass_frac_fuel > 1.0) ?
             1.0 :
             prob_parm.mass_frac_fuel);

      amrex::Real nonfuel_mass_fraction = 0.0;

      for (int n = 0; n < NUM_SPECIES; n++) {
        if (n != fuel_species_id) {
          nonfuel_mass_fraction += massfrac[n];
        }
      }

      for (int n = 0; n < NUM_SPECIES; n++) {
        if (n == fuel_species_id) {
          massfrac[n] = fuel_mass_fraction;
        } else if (nonfuel_mass_fraction > 0.0) {
          massfrac[n] *=
            (1.0 - fuel_mass_fraction) / nonfuel_mass_fraction;
        } else {
          massfrac[n] = 0.0;
        }
      }

      amrex::Real p_in =
        prob_parm.p_fuel_inlet;

      if (prob_parm.fuel_inlet_type == 1) {
        p_in =
          p_inter;
      }

      const amrex::Real T_in =
        prob_parm.T_fuel_inlet;

      eos.PYT2RE(
        p_in,
        massfrac,
        T_in,
        rho,
        eint);

      amrex::Real Cs_in;

      eos.RTY2Cs(
        rho,
        T_in,
        massfrac,
        Cs_in);

      if (prob_parm.fuel_inlet_type == 0) {
        u =
          -prob_parm.M_fuel_inlet * Cs_in;
      } else if (prob_parm.fuel_inlet_type == 1) {
        const amrex::Real pi =
          4.0 * std::atan(1.0);
        const amrex::Real fuel_port_area =
          pi *
          prob_parm.r_fuel_line_inner *
          prob_parm.r_fuel_line_inner;
        const amrex::Real design_inflow_speed =
          prob_parm.mdot_fuel_inlet / (rho * fuel_port_area);
        const amrex::Real backflow_guard_speed =
          (u > 0.0) ? u : 0.0;
        u =
          -((design_inflow_speed > backflow_guard_speed) ?
              design_inflow_speed :
              backflow_guard_speed);
      } else {
        amrex::Abort(
          "Invalid fuel inlet boundary type specified");
      }

      v = 0.0;
      w = 0.0;
      T = T_in;

    } else if (prob_parm.outlet_type == 0) {

      // First-order extrapolation.

    } else if (prob_parm.outlet_type == 1) {


      // One characteristic pressure-outlet state for every ghost cell.
      //
      // There is no x-dependent linear extrapolation and no switch to a
      // second thermodynamic state during backflow. The characteristic state
      // is constructed once, and its normal velocity is limited so that the
      // outlet cannot admit flow into the domain.

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

      amrex::Real cs_inter = 0.0;

      eos.RTY2Cs(
        rho_inter,
        T_inter,
        massfrac,
        cs_inter);

      const amrex::Real pb =
        prob_parm.p_exit;

      const amrex::Real ub =
        u_inter +
        (p_inter - pb) /
        (rho_inter * cs_inter);

      // Construct one EOS-consistent pressure-outlet state at p_exit while
      // retaining the interior temperature and composition.
      //
      // This same state is used for every outlet ghost cell, so there is no
      // linear extrapolation with ghost-cell depth.
      eos.PYT2RE(
        pb,
        massfrac,
        T_inter,
        rho,
        eint);

      // At the x-high boundary, positive u is outward.
      //
      // If the adjacent interior state points inward, require the exterior
      // velocity to be at least -u_inter. For an approximately symmetric
      // Riemann state, this gives a zero face-normal velocity.
      //
      // Otherwise, require the exterior velocity to be nonnegative.
      const amrex::Real u_min_no_backflow =
        (u_inter < 0.0) ?
          -u_inter :
          0.0;

      u =
        (ub > u_min_no_backflow) ?
          ub :
          u_min_no_backflow;

      v = v_inter;
      w = w_inter;
      T = T_inter;
      
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
