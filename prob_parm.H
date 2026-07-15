#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

AMREX_ENUM(
  TurbInflowType,
  None,
  TurbulentFluctuations,
  TotalVelocity);

struct ProbParmDevice
{
  // ---------------------------------------------------------------------------
  // Flow and thermodynamic parameters
  // ---------------------------------------------------------------------------

  amrex::Real M_inlet = 0.1;
  amrex::Real T_inlet = 300.0;
  amrex::Real p_inlet = 1.013e6;
  amrex::Real p_exit = 1.013e6;

  // Initial domain state.
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;

  // ---------------------------------------------------------------------------
  // Overall annular-sector envelope
  // ---------------------------------------------------------------------------

  amrex::Real L = 9.0;

  amrex::Real r_upper = 6.0;
  amrex::Real r_lower = 1.5;

  amrex::Real theta_deg = 20.0;

  // The x-low boundary is open only for r > inlet_lower_r.
  amrex::Real inlet_lower_r = 5.4;

  // Physical axial thickness of the inlet wall.
  amrex::Real inlet_cap_length = 0.25;

  // Nominal radial opening at the x-high boundary.
  amrex::Real outlet_upper_r = 4.0;
  amrex::Real outlet_lower_r = 2.6;

  // ---------------------------------------------------------------------------
  // U-shaped combustor liner
  // ---------------------------------------------------------------------------

  // Upstream dome axial interval.
  amrex::Real x_dome_lo = 1.40;
  amrex::Real x_dome_hi = 1.60;

  // Inner liner:
  //
  //   r_inner_liner_lo = inner-plenum side
  //   r_inner_liner_hi = combustion-chamber side
  //
  amrex::Real r_inner_liner_lo = 2.40;
  amrex::Real r_inner_liner_hi = 2.60;

  // Outer liner:
  //
  //   r_outer_liner_lo = combustion-chamber side
  //   r_outer_liner_hi = outer-plenum side
  //
  amrex::Real r_outer_liner_lo = 4.20;
  amrex::Real r_outer_liner_hi = 4.40;

  // Downstream radial turn of the outer liner.
  amrex::Real x_outer_turn_lo = 8.25;
  amrex::Real x_outer_turn_hi = 8.50;

  // Upper outlet lip radial interval.
  amrex::Real r_upper_exit_lo = 4.00;
  amrex::Real r_upper_exit_hi = 4.20;

  // ---------------------------------------------------------------------------
  // Inner casing walls
  // ---------------------------------------------------------------------------

  // Upstream inner-wall block.
  amrex::Real x_inner_wall_one_end = 1.00;
  amrex::Real r_inner_wall_one_top = 1.85;

  // Downstream inner-wall block.
  amrex::Real x_inner_wall_two_start = 5.50;
  amrex::Real r_inner_wall_two_top = 1.95;

  // Axial start of the outlet closures.
  //
  // The lower closure joins inner wall two to the inner liner.
  // The upper closure joins the upper outlet lip to the outer casing.
  //
  amrex::Real x_outlet_wall_lo = 8.70;

  // ---------------------------------------------------------------------------
  // Boundary and sponge options
  // ---------------------------------------------------------------------------

  int inlet_type = 0;
  int outlet_type = 1;

  bool do_sponge_zones = true;

  TurbInflowType turb_inflow_type =
    TurbInflowType::None;
};

struct ProbParmHost
{
  ProbParmHost() = default;
};

#endif
