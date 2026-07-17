#include "prob.H"

#include <AMReX_MFIter.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Vector.H>

#include <cmath>

namespace
{
void
validate_geometry(const ProbParmDevice& pp)
{
  auto require = [](const bool condition, const char* message) {
    if (!condition) {
      amrex::Abort(message);
    }
  };

  require(pp.L > 0.0, "prob.L must be positive");
  require(
    pp.r_lower > 0.0 && pp.r_lower < pp.r_upper,
    "Require 0 < prob.r_lower < prob.r_upper");
  require(
    pp.theta_deg > 0.0 && pp.theta_deg < 180.0,
    "Require 0 < prob.theta_deg < 180");

  require(
    pp.inlet_cap_length > 0.0 && pp.inlet_cap_length < pp.x_dome_lo,
    "Require 0 < prob.inlet_cap_length < prob.x_dome_lo");
  require(
    pp.inlet_lower_r > pp.r_lower && pp.inlet_lower_r < pp.r_upper,
    "prob.inlet_lower_r must lie inside the annulus");
  require(
    pp.fuel_species_id >= -1 && pp.fuel_species_id < NUM_SPECIES,
    "prob.fuel_species_id must be -1 for mechanism default or name a valid "
    "species");
  require(
    pp.mass_frac_fuel >= 0.0 && pp.mass_frac_fuel <= 1.0,
    "Require 0 <= prob.mass_frac_fuel <= 1");
  require(
    pp.air_mole_frac_o2 >= 0.0 && pp.air_mole_frac_n2 >= 0.0 &&
      pp.air_mole_frac_o2 + pp.air_mole_frac_n2 > 0.0,
    "Air O2/N2 mole fractions must be nonnegative with positive sum");
  require(
    pp.secondary_inlet_split_r > pp.r_lower &&
      pp.secondary_inlet_split_r < pp.r_upper,
    "prob.secondary_inlet_split_r must lie inside the annulus");

  require(
    0.0 <= pp.x_inlet_vane_lo &&
      pp.x_inlet_vane_lo < pp.x_inlet_vane_hi &&
      pp.x_inlet_vane_hi < pp.x_dome_lo,
    "Require 0 <= x_inlet_vane_lo < x_inlet_vane_hi < x_dome_lo");
  require(
    pp.r_lower < pp.r_inlet_vane_lo &&
      pp.r_inlet_vane_lo < pp.r_inlet_vane_hi &&
      pp.r_inlet_vane_hi <= pp.inlet_lower_r,
    "Inlet vane radii must be ordered below or at prob.inlet_lower_r");

  require(
    0.0 <= pp.x_evaporator_lo &&
      pp.x_evaporator_lo < pp.x_evaporator_bore_hi &&
      pp.x_evaporator_bore_hi <= pp.x_outer_turn_hi &&
      pp.x_outer_turn_hi < pp.L,
    "Require ordered evaporator x extents inside the domain");
  require(
    pp.r_evaporator_bore > 0.0 &&
      pp.r_evaporator_bore < pp.r_evaporator_outer,
    "Require 0 < r_evaporator_bore < r_evaporator_outer");

  require(
    pp.x_fuel_line_lo >= pp.x_evaporator_lo &&
      pp.x_fuel_line_lo < pp.L,
    "Require x_evaporator_lo <= x_fuel_line_lo < L");
  require(
    pp.r_fuel_line_inner > 0.0 &&
      pp.r_fuel_line_inner < pp.r_fuel_line_outer &&
      pp.r_fuel_line_outer < pp.r_evaporator_bore,
    "Require 0 < r_fuel_line_inner < r_fuel_line_outer < r_evaporator_bore");

  if (pp.do_spark) {
    require(
      pp.x_spark > 0.0 && pp.x_spark < pp.L,
      "prob.x_spark must lie inside the domain when prob.do_spark = true");
    require(
      pp.r_spark > pp.r_lower && pp.r_spark < pp.r_upper,
      "prob.r_spark must lie inside the annulus when prob.do_spark = true");
    require(
      pp.spark_radius > 0.0,
      "prob.spark_radius must be positive when prob.do_spark = true");
    require(
      pp.T_spark > 0.0,
      "prob.T_spark must be positive when prob.do_spark = true");
  }

  for (int n = 0; n < ProbParmDevice::num_liner_hole_rows; ++n) {
    require(
      pp.liner_hole_x[n] > 0.0 && pp.liner_hole_x[n] < pp.L,
      "Liner hole x-stations must lie inside the domain");
    require(
      pp.liner_hole_radius[n] > 0.0,
      "Liner hole radii must be positive");
    require(
      pp.liner_hole_spacing_deg[n] > 0.0,
      "Liner hole spacing must be positive");
    require(
      pp.liner_hole_theta_start_deg[n] >= 0.0 &&
        pp.liner_hole_theta_start_deg[n] <= pp.theta_deg,
      "Liner hole theta_start must lie inside the sector");
    require(
      pp.liner_hole_omit_every_n[n] >= 0,
      "Liner hole omit cadence must be nonnegative");
  }

  require(
    0.0 < pp.x_dome_lo && pp.x_dome_lo < pp.x_dome_hi &&
      pp.x_dome_hi < pp.L,
    "Require 0 < x_dome_lo < x_dome_hi < L");

  require(
    pp.r_lower < pp.r_inner_liner_lo &&
      pp.r_inner_liner_lo < pp.r_inner_liner_hi &&
      pp.r_inner_liner_hi < pp.r_outer_liner_lo &&
      pp.r_outer_liner_lo < pp.r_outer_liner_hi &&
      pp.r_outer_liner_hi < pp.r_upper,
    "Liner radii must be ordered from the inner casing to the outer casing");

  require(
    pp.x_dome_hi < pp.x_outer_turn_lo &&
      pp.x_outer_turn_lo < pp.x_outer_turn_hi &&
      pp.x_outer_turn_hi < pp.L,
    "Require x_dome_hi < x_outer_turn_lo < x_outer_turn_hi < L");

  require(
    pp.r_inner_liner_hi < pp.r_upper_exit_lo &&
      pp.r_upper_exit_lo < pp.r_upper_exit_hi &&
      pp.r_upper_exit_hi <= pp.r_outer_liner_lo,
    "Upper exit lip must lie between the chamber and the outer liner");

  require(
    0.0 < pp.x_inner_wall_one_end &&
      pp.x_inner_wall_one_end < pp.x_dome_lo,
    "Inner wall one must end upstream of the dome");
  require(
    pp.r_lower < pp.r_inner_wall_one_top &&
      pp.r_inner_wall_one_top < pp.r_inner_liner_lo,
    "Inner wall one must remain below the inner liner");

  require(
    pp.x_dome_hi < pp.x_inner_wall_two_start &&
      pp.x_inner_wall_two_start < pp.x_outlet_wall_lo &&
      pp.x_outlet_wall_lo < pp.L,
    "Require x_dome_hi < x_inner_wall_two_start < x_outlet_wall_lo < L");
  require(
    pp.r_lower < pp.r_inner_wall_two_top &&
      pp.r_inner_wall_two_top < pp.r_inner_liner_lo,
    "Inner wall two must remain below the inner liner");

  require(
    pp.outlet_lower_r >= pp.r_inner_liner_hi &&
      pp.outlet_upper_r <= pp.r_upper_exit_lo &&
      pp.outlet_lower_r < pp.outlet_upper_r,
    "Nominal outlet radii must lie inside the geometric outlet opening");
}
} // namespace

void
pc_prob_close()
{
}

void
parse_params(ProbParmDevice* prob_parm_device)
{
  amrex::ParmParse pp("prob");

  pp.query("M_inlet", prob_parm_device->M_inlet);
  pp.query("T_inlet", prob_parm_device->T_inlet);
  pp.query("p_inlet", prob_parm_device->p_inlet);
  pp.query("p_exit", prob_parm_device->p_exit);
  pp.query("M_fuel_inlet", prob_parm_device->M_fuel_inlet);
  pp.query("T_fuel_inlet", prob_parm_device->T_fuel_inlet);
  pp.query("p_fuel_inlet", prob_parm_device->p_fuel_inlet);
  pp.query("fuel_species_id", prob_parm_device->fuel_species_id);
  pp.query("mass_frac_fuel", prob_parm_device->mass_frac_fuel);
  pp.query("debug_fuel_port", prob_parm_device->debug_fuel_port);
  pp.query("air_mole_frac_o2", prob_parm_device->air_mole_frac_o2);
  pp.query("air_mole_frac_n2", prob_parm_device->air_mole_frac_n2);
  pp.query("p0", prob_parm_device->p0);
  pp.query("T0", prob_parm_device->T0);
  pp.query("inlet_type", prob_parm_device->inlet_type);
  pp.query("outlet_type", prob_parm_device->outlet_type);
  pp.query("do_sponge_zones", prob_parm_device->do_sponge_zones);
  pp.query("turb_inflow_type", prob_parm_device->turb_inflow_type);

  // Overall annular-sector geometry.
  pp.query("L", prob_parm_device->L);
  pp.query("r_upper", prob_parm_device->r_upper);
  pp.query("r_lower", prob_parm_device->r_lower);
  pp.query("theta_deg", prob_parm_device->theta_deg);

  pp.query("inlet_lower_r", prob_parm_device->inlet_lower_r);
  pp.query("inlet_cap_length", prob_parm_device->inlet_cap_length);
  pp.query(
    "secondary_inlet_split_r",
    prob_parm_device->secondary_inlet_split_r);
  pp.query("outlet_upper_r", prob_parm_device->outlet_upper_r);
  pp.query("outlet_lower_r", prob_parm_device->outlet_lower_r);

  // Inlet vane.
  pp.query("x_inlet_vane_lo", prob_parm_device->x_inlet_vane_lo);
  pp.query("x_inlet_vane_hi", prob_parm_device->x_inlet_vane_hi);
  pp.query("r_inlet_vane_lo", prob_parm_device->r_inlet_vane_lo);
  pp.query("r_inlet_vane_hi", prob_parm_device->r_inlet_vane_hi);

  // Evaporator.
  pp.query("x_evaporator_lo", prob_parm_device->x_evaporator_lo);
  pp.query("x_evaporator_bore_hi", prob_parm_device->x_evaporator_bore_hi);
  pp.query(
    "y_evaporator_center",
    prob_parm_device->y_evaporator_center);
  pp.query(
    "z_evaporator_center",
    prob_parm_device->z_evaporator_center);
  pp.query("r_evaporator_outer", prob_parm_device->r_evaporator_outer);
  pp.query("r_evaporator_bore", prob_parm_device->r_evaporator_bore);

  // Fuel line.
  pp.query("x_fuel_line_lo", prob_parm_device->x_fuel_line_lo);
  pp.query("r_fuel_line_outer", prob_parm_device->r_fuel_line_outer);
  pp.query("r_fuel_line_inner", prob_parm_device->r_fuel_line_inner);

  // Spark ignition. The spark center is specified in (x, cylindrical-r), and
  // spark_radius controls the meridional radius of the heated restart region.
  pp.query("do_spark", prob_parm_device->do_spark);
  pp.query("x_spark", prob_parm_device->x_spark);
  pp.query("r_spark", prob_parm_device->r_spark);
  pp.query("spark_radius", prob_parm_device->spark_radius);
  pp.query("T_spark", prob_parm_device->T_spark);

  auto query_real_array =
    [&](const char* name, amrex::Real* values) {
      amrex::Vector<amrex::Real> queried_values;
      if (
        pp.queryarr(
          name,
          queried_values,
          0,
          ProbParmDevice::num_liner_hole_rows)) {
        for (int n = 0; n < ProbParmDevice::num_liner_hole_rows; ++n) {
          values[n] = queried_values[n];
        }
      }
    };

  auto query_int_array =
    [&](const char* name, int* values) {
      amrex::Vector<int> queried_values;
      if (
        pp.queryarr(
          name,
          queried_values,
          0,
          ProbParmDevice::num_liner_hole_rows)) {
        for (int n = 0; n < ProbParmDevice::num_liner_hole_rows; ++n) {
          values[n] = queried_values[n];
        }
      }
    };

  // Liner holes.
  query_real_array("liner_hole_x", prob_parm_device->liner_hole_x);
  query_real_array("liner_hole_radius", prob_parm_device->liner_hole_radius);
  query_real_array(
    "liner_hole_spacing_deg",
    prob_parm_device->liner_hole_spacing_deg);
  query_real_array(
    "liner_hole_theta_start_deg",
    prob_parm_device->liner_hole_theta_start_deg);
  query_int_array(
    "liner_hole_omit_every_n",
    prob_parm_device->liner_hole_omit_every_n);

  // U-shaped combustor liner.
  pp.query("x_dome_lo", prob_parm_device->x_dome_lo);
  pp.query("x_dome_hi", prob_parm_device->x_dome_hi);

  pp.query("r_inner_liner_lo", prob_parm_device->r_inner_liner_lo);
  pp.query("r_inner_liner_hi", prob_parm_device->r_inner_liner_hi);
  pp.query("r_outer_liner_lo", prob_parm_device->r_outer_liner_lo);
  pp.query("r_outer_liner_hi", prob_parm_device->r_outer_liner_hi);

  pp.query("x_outer_turn_lo", prob_parm_device->x_outer_turn_lo);
  pp.query("x_outer_turn_hi", prob_parm_device->x_outer_turn_hi);
  pp.query("r_upper_exit_lo", prob_parm_device->r_upper_exit_lo);
  pp.query("r_upper_exit_hi", prob_parm_device->r_upper_exit_hi);

  pp.query(
    "x_inner_wall_one_end", prob_parm_device->x_inner_wall_one_end);
  pp.query(
    "r_inner_wall_one_top", prob_parm_device->r_inner_wall_one_top);

  pp.query(
    "x_inner_wall_two_start", prob_parm_device->x_inner_wall_two_start);
  pp.query(
    "r_inner_wall_two_top", prob_parm_device->r_inner_wall_two_top);

  pp.query("x_outlet_wall_lo", prob_parm_device->x_outlet_wall_lo);

  validate_geometry(*prob_parm_device);
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  parse_params(PeleC::h_prob_parm_device);
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
  parse_params(PeleC::h_prob_parm_device);

  const ProbParmDevice pp =
    *PeleC::h_prob_parm_device;

  if (!pp.do_spark) {
    return;
  }

  auto& state =
    get_new_data(State_Type);

  const amrex::GeometryData geomdata =
    Geom().data();

  const amrex::Real x_spark =
    pp.x_spark;
  const amrex::Real r_spark =
    pp.r_spark;
  const amrex::Real spark_radius =
    pp.spark_radius;
  const amrex::Real T_spark =
    pp.T_spark;

  for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx =
      mfi.tilebox();

    auto const state_arr =
      state.array(mfi);

    amrex::ParallelFor(
      bx,
      [=] AMREX_GPU_DEVICE(const int i, const int j, const int k) noexcept {
        const amrex::Real* prob_lo =
          geomdata.ProbLo();
        const amrex::Real* dx =
          geomdata.CellSize();

        const amrex::Real x =
          prob_lo[0] + (static_cast<amrex::Real>(i) + 0.5) * dx[0];
        const amrex::Real y =
          prob_lo[1] + (static_cast<amrex::Real>(j) + 0.5) * dx[1];
        const amrex::Real z =
          prob_lo[2] + (static_cast<amrex::Real>(k) + 0.5) * dx[2];

        const amrex::Real r =
          std::sqrt(y * y + z * z);

        const amrex::Real dx_spark =
          x - x_spark;
        const amrex::Real dr_spark =
          r - r_spark;

        if (
          dx_spark * dx_spark + dr_spark * dr_spark <=
          spark_radius * spark_radius) {

          const amrex::Real rho =
            state_arr(i, j, k, URHO);

          amrex::Real massfrac[NUM_SPECIES];

          for (int n = 0; n < NUM_SPECIES; n++) {
            massfrac[n] =
              state_arr(i, j, k, UFS + n) / rho;
          }

          auto eos =
            pele::physics::PhysicsType::eos();

          amrex::Real eint = 0.0;

          eos.RTY2E(
            rho,
            T_spark,
            massfrac,
            eint);

          const amrex::Real u =
            state_arr(i, j, k, UMX) / rho;
          const amrex::Real v =
            state_arr(i, j, k, UMY) / rho;
          const amrex::Real w =
            state_arr(i, j, k, UMZ) / rho;

          state_arr(i, j, k, UEINT) =
            rho * eint;
          state_arr(i, j, k, UEDEN) =
            rho *
            (
              eint +
              0.5 * (u * u + v * v + w * w)
            );
          state_arr(i, j, k, UTEMP) =
            T_spark;
        }
      });
  }
}

void
EBAnnularSector::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  parse_params(PeleC::h_prob_parm_device);
  ProbParmDevice const* pp = PeleC::h_prob_parm_device;

  const amrex::Real L = pp->L;
  const amrex::Real theta =
    pp->theta_deg * static_cast<amrex::Real>(M_PI) / 180.0;

  // All cylindrical surfaces are coaxial with the x-axis. The x-coordinate
  // of this point does not affect an axis-0 cylinder, but retaining 0.5*L
  // follows the convention used by the original problem.
  const amrex::RealArray axis_point{
    AMREX_D_DECL(0.5 * L, 0.0, 0.0)};

  auto make_capped_cylinder =
    [&](const amrex::Real xlo,
        const amrex::Real xhi,
        const amrex::Real radius,
        const amrex::Real y_center,
        const amrex::Real z_center) {
      const amrex::RealArray cylinder_axis_point{
        AMREX_D_DECL(0.5 * (xlo + xhi), y_center, z_center)};

      amrex::EB2::CylinderIF inside_radius(
        radius, 0, cylinder_axis_point, false);

      amrex::EB2::PlaneIF right_of_xlo(
        {AMREX_D_DECL(xlo, 0.0, 0.0)},
        {AMREX_D_DECL(1.0, 0.0, 0.0)},
        true);

      amrex::EB2::PlaneIF left_of_xhi(
        {AMREX_D_DECL(xhi, 0.0, 0.0)},
        {AMREX_D_DECL(-1.0, 0.0, 0.0)},
        true);

      auto axial_slab =
        amrex::EB2::makeIntersection(
          right_of_xlo,
          left_of_xhi);

      return amrex::EB2::makeIntersection(
        inside_radius,
        axial_slab);
    };

  auto make_liner_hole =
    [&](const amrex::Real x,
        const amrex::Real radius,
        const amrex::Real liner_radius_mid,
        const amrex::Real liner_thickness,
        const amrex::Real theta_hole) {
      const amrex::Real height =
        liner_thickness + 2.0 * radius;

      const amrex::RealArray hole_center{
        AMREX_D_DECL(x, liner_radius_mid, 0.0)};

      amrex::EB2::CylinderIF hole(
        radius, height, 1, hole_center, false);

      return amrex::EB2::rotate(hole, theta_hole, 0);
    };

  // Construct a solid annular block:
  //
  //   xlo < x < xhi,
  //   rlo < sqrt(y^2 + z^2) < rhi.
  //
  auto make_annular_block =
    [&](const amrex::Real xlo,
        const amrex::Real xhi,
        const amrex::Real rlo,
        const amrex::Real rhi) {
      // Solid for r > rlo.
      amrex::EB2::CylinderIF outside_rlo(
        rlo, 0, axis_point, true);

      // Solid for r < rhi.
      amrex::EB2::CylinderIF inside_rhi(
        rhi, 0, axis_point, false);

      auto radial_shell =
        amrex::EB2::makeIntersection(
          outside_rlo,
          inside_rhi);

      // Solid for x > xlo.
      amrex::EB2::PlaneIF right_of_xlo(
        {AMREX_D_DECL(xlo, 0.0, 0.0)},
        {AMREX_D_DECL(1.0, 0.0, 0.0)},
        true);

      // Solid for x < xhi.
      amrex::EB2::PlaneIF left_of_xhi(
        {AMREX_D_DECL(xhi, 0.0, 0.0)},
        {AMREX_D_DECL(-1.0, 0.0, 0.0)},
        true);

      auto axial_slab =
        amrex::EB2::makeIntersection(
          right_of_xlo,
          left_of_xhi);

      return amrex::EB2::makeIntersection(
        radial_shell,
        axial_slab);
    };

  // ---------------------------------------------------------------------------
  // Original annular-sector envelope
  // ---------------------------------------------------------------------------

  // r < r_lower is solid.
  amrex::EB2::CylinderIF inner_solid(
    pp->r_lower,
    0,
    axis_point,
    false);

  // r > r_upper is solid.
  amrex::EB2::CylinderIF outer_solid(
    pp->r_upper,
    0,
    axis_point,
    true);

  // High-theta side is solid. The other angular side is expected to coincide
  // with a Cartesian domain boundary.
  const amrex::Real c = std::cos(theta);
  const amrex::Real s = std::sin(theta);

  amrex::EB2::PlaneIF sector_solid(
    {AMREX_D_DECL(0.0, 0.0, 0.0)},
    {AMREX_D_DECL(0.0, -s, c)},
    true);

  // The x-low boundary is open only above inlet_lower_r.
  auto inlet_wall_block = make_annular_block(
    0.0,
    pp->inlet_cap_length,
    pp->r_lower,
    pp->inlet_lower_r);

  auto inlet_vane = make_annular_block(
    pp->x_inlet_vane_lo,
    pp->x_inlet_vane_hi,
    pp->r_inlet_vane_lo,
    pp->r_inlet_vane_hi);

  auto evaporator_outer = make_capped_cylinder(
    pp->x_evaporator_lo,
    pp->x_outer_turn_hi,
    pp->r_evaporator_outer,
    pp->y_evaporator_center,
    pp->z_evaporator_center);

  auto evaporator_bore = make_capped_cylinder(
    pp->x_evaporator_lo,
    pp->x_evaporator_bore_hi,
    pp->r_evaporator_bore,
    pp->y_evaporator_center,
    pp->z_evaporator_center);

  auto evaporator_bore_cutout =
    amrex::EB2::makeComplement(evaporator_bore);

  auto fuel_line_outer = make_capped_cylinder(
    pp->x_fuel_line_lo,
    L,
    pp->r_fuel_line_outer,
    pp->y_evaporator_center,
    pp->z_evaporator_center);

  auto fuel_line_bore = make_capped_cylinder(
    pp->x_fuel_line_lo,
    L,
    pp->r_fuel_line_inner,
    pp->y_evaporator_center,
    pp->z_evaporator_center);

  const amrex::Real outer_liner_mid =
    0.5 * (pp->r_outer_liner_lo + pp->r_outer_liner_hi);
  const amrex::Real outer_liner_thickness =
    pp->r_outer_liner_hi - pp->r_outer_liner_lo;

  const amrex::Real inner_liner_mid =
    0.5 * (pp->r_inner_liner_lo + pp->r_inner_liner_hi);
  const amrex::Real inner_liner_thickness =
    pp->r_inner_liner_hi - pp->r_inner_liner_lo;

  constexpr amrex::Real deg_to_rad =
    static_cast<amrex::Real>(M_PI) / 180.0;

  auto outer_hole_1_0 = make_liner_hole(
    pp->liner_hole_x[0], pp->liner_hole_radius[0], outer_liner_mid,
    outer_liner_thickness, pp->liner_hole_theta_start_deg[0] * deg_to_rad);
  auto outer_hole_1_1 = make_liner_hole(
    pp->liner_hole_x[0], pp->liner_hole_radius[0], outer_liner_mid,
    outer_liner_thickness,
    (
      pp->liner_hole_theta_start_deg[0] +
      pp->liner_hole_spacing_deg[0]
    ) * deg_to_rad);
  auto outer_hole_1_2 = make_liner_hole(
    pp->liner_hole_x[0], pp->liner_hole_radius[0], outer_liner_mid,
    outer_liner_thickness,
    (
      pp->liner_hole_theta_start_deg[0] +
      2.0 * pp->liner_hole_spacing_deg[0]
    ) * deg_to_rad);
  auto outer_hole_2 = make_liner_hole(
    pp->liner_hole_x[1], pp->liner_hole_radius[1], outer_liner_mid,
    outer_liner_thickness, pp->liner_hole_theta_start_deg[1] * deg_to_rad);
  auto outer_hole_3 = make_liner_hole(
    pp->liner_hole_x[2], pp->liner_hole_radius[2], outer_liner_mid,
    outer_liner_thickness, pp->liner_hole_theta_start_deg[2] * deg_to_rad);
  auto outer_hole_4 = make_liner_hole(
    pp->liner_hole_x[3], pp->liner_hole_radius[3], outer_liner_mid,
    outer_liner_thickness, pp->liner_hole_theta_start_deg[3] * deg_to_rad);
  auto outer_hole_5 = make_liner_hole(
    pp->liner_hole_x[4], pp->liner_hole_radius[4], outer_liner_mid,
    outer_liner_thickness, pp->liner_hole_theta_start_deg[4] * deg_to_rad);
  auto outer_hole_6 = make_liner_hole(
    pp->liner_hole_x[5], pp->liner_hole_radius[5], outer_liner_mid,
    outer_liner_thickness,
    (
      pp->liner_hole_theta_start_deg[5] +
      pp->liner_hole_spacing_deg[5]
    ) * deg_to_rad);
  auto outer_hole_7 = make_liner_hole(
    pp->liner_hole_x[6], pp->liner_hole_radius[6], outer_liner_mid,
    outer_liner_thickness, pp->liner_hole_theta_start_deg[6] * deg_to_rad);

  auto inner_hole_1_0 = make_liner_hole(
    pp->liner_hole_x[0], pp->liner_hole_radius[0], inner_liner_mid,
    inner_liner_thickness, pp->liner_hole_theta_start_deg[0] * deg_to_rad);
  auto inner_hole_1_1 = make_liner_hole(
    pp->liner_hole_x[0], pp->liner_hole_radius[0], inner_liner_mid,
    inner_liner_thickness,
    (
      pp->liner_hole_theta_start_deg[0] +
      pp->liner_hole_spacing_deg[0]
    ) * deg_to_rad);
  auto inner_hole_1_2 = make_liner_hole(
    pp->liner_hole_x[0], pp->liner_hole_radius[0], inner_liner_mid,
    inner_liner_thickness,
    (
      pp->liner_hole_theta_start_deg[0] +
      2.0 * pp->liner_hole_spacing_deg[0]
    ) * deg_to_rad);
  auto inner_hole_2 = make_liner_hole(
    pp->liner_hole_x[1], pp->liner_hole_radius[1], inner_liner_mid,
    inner_liner_thickness, pp->liner_hole_theta_start_deg[1] * deg_to_rad);
  auto inner_hole_3 = make_liner_hole(
    pp->liner_hole_x[2], pp->liner_hole_radius[2], inner_liner_mid,
    inner_liner_thickness, pp->liner_hole_theta_start_deg[2] * deg_to_rad);
  auto inner_hole_4 = make_liner_hole(
    pp->liner_hole_x[3], pp->liner_hole_radius[3], inner_liner_mid,
    inner_liner_thickness, pp->liner_hole_theta_start_deg[3] * deg_to_rad);
  auto inner_hole_5 = make_liner_hole(
    pp->liner_hole_x[4], pp->liner_hole_radius[4], inner_liner_mid,
    inner_liner_thickness, pp->liner_hole_theta_start_deg[4] * deg_to_rad);
  auto inner_hole_6 = make_liner_hole(
    pp->liner_hole_x[5], pp->liner_hole_radius[5], inner_liner_mid,
    inner_liner_thickness,
    (
      pp->liner_hole_theta_start_deg[5] +
      pp->liner_hole_spacing_deg[5]
    ) * deg_to_rad);
  auto inner_hole_7 = make_liner_hole(
    pp->liner_hole_x[6], pp->liner_hole_radius[6], inner_liner_mid,
    inner_liner_thickness, pp->liner_hole_theta_start_deg[6] * deg_to_rad);

  auto liner_holes =
    amrex::EB2::makeUnion(
      outer_hole_1_0,
      outer_hole_1_1,
      outer_hole_1_2,
      outer_hole_2,
      outer_hole_3,
      outer_hole_4,
      outer_hole_5,
      outer_hole_6,
      outer_hole_7,
      inner_hole_1_0,
      inner_hole_1_1,
      inner_hole_1_2,
      inner_hole_2,
      inner_hole_3,
      inner_hole_4,
      inner_hole_5,
      inner_hole_6,
      inner_hole_7);

  // ---------------------------------------------------------------------------
  // Red combustor walls from the meridional sketch
  // ---------------------------------------------------------------------------

  // Points 1-2-8-9: outer cylindrical liner.
  auto outer_liner = make_annular_block(
    pp->x_dome_lo,
    pp->x_outer_turn_hi,
    pp->r_outer_liner_lo,
    pp->r_outer_liner_hi);

  // Points 1-9-10-13: upstream annular dome joining both liners.
  auto upstream_dome = make_annular_block(
    pp->x_dome_lo,
    pp->x_dome_hi,
    pp->r_inner_liner_lo,
    pp->r_outer_liner_hi);

  // Points 10-11-12-13: inner cylindrical liner.
  auto inner_liner = make_annular_block(
    pp->x_dome_lo,
    L,
    pp->r_inner_liner_lo,
    pp->r_inner_liner_hi);

  // Points 2-3-7-8: downstream radial turn of the outer liner.
  auto outer_exit_turn = make_annular_block(
    pp->x_outer_turn_lo,
    pp->x_outer_turn_hi,
    pp->r_upper_exit_lo,
    pp->r_outer_liner_hi);

  // Points 3-4-5-6: upper axial outlet lip. Keep this wall upstream of
  // x_outlet_wall_lo so r > secondary_inlet_split_r reaches the x-high face.
  auto upper_exit_lip = make_annular_block(
    pp->x_outer_turn_lo,
    pp->x_outlet_wall_lo,
    pp->r_upper_exit_lo,
    pp->r_upper_exit_hi);

  // Points 14-15-16-17: inner wall one.
  auto inner_wall_one = make_annular_block(
    0.0,
    pp->x_inner_wall_one_end,
    pp->r_lower,
    pp->r_inner_wall_one_top);

  // Points 18-19-22-23: inner wall two.
  auto inner_wall_two = make_annular_block(
    pp->x_inner_wall_two_start,
    L,
    pp->r_lower,
    pp->r_inner_wall_two_top);

  // Close the inner plenum at the outlet.
  auto lower_outlet_wall = make_annular_block(
    pp->x_outlet_wall_lo,
    L,
    pp->r_inner_wall_two_top,
    pp->r_inner_liner_lo);

  // Close the upper x-high face. The fuel-line bore cutout below is applied
  // later, leaving a small fuel-inlet port behind the evaporator.
  auto upper_fuel_inlet_wall = make_annular_block(
    pp->x_outlet_wall_lo,
    L,
    pp->secondary_inlet_split_r,
    pp->r_upper);

  auto combustor_walls_before_evaporator_bore =
    amrex::EB2::makeUnion(
      outer_liner,
      upstream_dome,
      inner_liner,
      outer_exit_turn,
      upper_exit_lip,
      inlet_vane,
      evaporator_outer,
      inner_wall_one,
      inner_wall_two,
      lower_outlet_wall,
      upper_fuel_inlet_wall);

  auto solid_before_evaporator_bore =
    amrex::EB2::makeUnion(
      inner_solid,
      outer_solid,
      inlet_wall_block,
      combustor_walls_before_evaporator_bore);

  // Apply the evaporator bore before adding the fuel-line wall so the
  // evaporator interior cuts through the liner but does not erase the
  // dedicated fuel-line tube.
  auto solid_after_evaporator_bore =
    amrex::EB2::makeIntersection(
      solid_before_evaporator_bore,
      evaporator_bore_cutout);

  auto solid_before_cutouts =
    amrex::EB2::makeUnion(
      solid_after_evaporator_bore,
      fuel_line_outer);

  auto fuel_line_bore_cutout =
    amrex::EB2::makeComplement(fuel_line_bore);

  auto liner_hole_cutouts =
    amrex::EB2::makeComplement(liner_holes);

  auto cut_solid_region =
    amrex::EB2::makeIntersection(
      solid_before_cutouts,
      fuel_line_bore_cutout,
      liner_hole_cutouts);

  // Apply the sector wall after subtracting holes so hole cutouts cannot
  // perforate the symmetry-sector boundary.
  auto solid_region =
    amrex::EB2::makeUnion(
      cut_solid_region,
      sector_solid);

  auto gshop =
    amrex::EB2::makeShop(solid_region);

  amrex::EB2::Build(
    gshop,
    geom,
    max_coarsening_level,
    max_coarsening_level,
    4,
    false);
}
