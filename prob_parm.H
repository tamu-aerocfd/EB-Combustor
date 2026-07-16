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

  amrex::Real M_fuel_inlet = 0.1;
  amrex::Real T_fuel_inlet = 300.0;
  amrex::Real p_fuel_inlet = 1.013e6;
  int fuel_species_id = 0;
  amrex::Real mass_frac_fuel = 1.0;

  // Initial domain state.
  amrex::Real p0 = 1.013e6; // [erg cm^-3]
  amrex::Real T0 = 300.0;

  // ---------------------------------------------------------------------------
  // Overall annular-sector envelope
  // ---------------------------------------------------------------------------

  amrex::Real L = 9.0;

  amrex::Real r_upper = 6.0;
  amrex::Real r_lower = 1.5;

  amrex::Real theta_deg = 15.0;

  // The x-low boundary is open only for r > inlet_lower_r.
  amrex::Real inlet_lower_r = 5.4;

  // Physical axial thickness of the inlet wall.
  amrex::Real inlet_cap_length = 0.25;

  // Split radius at the x-high boundary: r below this value is outlet,
  // r above this value is the secondary fuel inlet.
  amrex::Real secondary_inlet_split_r = 4.0;

  // Nominal radial opening at the x-high outlet.
  amrex::Real outlet_upper_r = 4.0;
  amrex::Real outlet_lower_r = 2.6;

  // ---------------------------------------------------------------------------
  // Inlet vane
  // ---------------------------------------------------------------------------

  amrex::Real x_inlet_vane_lo = 0.0;
  amrex::Real x_inlet_vane_hi = 0.5;
  amrex::Real r_inlet_vane_lo = 5.3;
  amrex::Real r_inlet_vane_hi = 5.4;

  // ---------------------------------------------------------------------------
  // Evaporator
  // ---------------------------------------------------------------------------

  amrex::Real x_evaporator_lo = 2.0;
  amrex::Real x_evaporator_bore_hi = 8.5;
  amrex::Real y_evaporator_center = 4.6;
  amrex::Real z_evaporator_center = 0.0;
  amrex::Real r_evaporator_outer = 0.4;
  amrex::Real r_evaporator_bore = 0.3;

  // ---------------------------------------------------------------------------
  // Fuel line
  // ---------------------------------------------------------------------------

  amrex::Real x_fuel_line_lo = 2.0;
  amrex::Real r_fuel_line_outer = 0.20;
  amrex::Real r_fuel_line_inner = 0.15;

  // ---------------------------------------------------------------------------
  // Spark ignition
  // ---------------------------------------------------------------------------

  bool do_spark = false;
  amrex::Real x_spark = 4.5;
  amrex::Real r_spark = 3.5;
  amrex::Real spark_radius = 0.1;
  amrex::Real T_spark = 2000.0;

  // ---------------------------------------------------------------------------
  // Liner holes
  // ---------------------------------------------------------------------------

  static constexpr int num_liner_hole_rows = 7;

  amrex::Real liner_hole_x[num_liner_hole_rows] =
    {2.3, 3.05, 3.9, 4.7, 5.6, 6.5, 7.5};
  amrex::Real liner_hole_radius[num_liner_hole_rows] =
    {0.1, 0.15, 0.2, 0.2, 0.2, 0.25, 0.1};
  amrex::Real liner_hole_spacing_deg[num_liner_hole_rows] =
    {6.0, 30.0, 30.0, 30.0, 30.0, 9.0, 15.0};
  amrex::Real liner_hole_theta_start_deg[num_liner_hole_rows] =
    {0.0, 15.0, 15.0, 15.0, 15.0, 0.0, 7.5};
  int liner_hole_omit_every_n[num_liner_hole_rows] =
    {0, 0, 0, 0, 0, 3, 0};

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
  // The upper side remains open at x-high for the secondary fuel inlet.
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
