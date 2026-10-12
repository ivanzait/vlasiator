#pragma once
/*
 * Adaptive orthonormal Hermite transform of a dense 3D velocity cube.
 *
 * Self-contained (no Vlasiator/dccrg/MPI types). Port of the Python reference
 * faiser/data_processing/vdf_tools.py (adaptive_transform & friends).
 * All computation is done in double precision.
 *
 * LAYOUT AND UNITS
 *   cube[(iz * ny + iy) * nx + ix]: C order, z slowest, x fastest (cube[iz][iy][ix]).
 *   Cell centres: v_d[i] = v0[d] + i * dv[d] (struct Grid). Velocities are in m/s, f in its own units
 *   (s^3/m^6 in Vlasiator). Coefficients are plain sums over cell centres,
 *   c = dv_x dv_y dv_z * sum f * phi_l(z) phi_m(y) phi_n(x), i.e. in units of f * (m/s)^(3/2). Hence a cube that is
 *   zero outside a sub-box may be cropped to that box without changing any result.
 *   The cubic API taking (vlim, vlen) is a thin wrapper over
 *   Grid{v0 = -vlim + dv/2, dv = 2*vlim/vlen, shape = vlen^3}.
 *
 * BASIS
 *   One isotropic thermal velocity vth (the standard deviation of the Gaussian) and the bulk velocity u
 *   define the 1-D orthonormal Hermite functions phi_n((v - u_d) / vth) / sqrt(vth), n = 0, 1, ..., used
 *   along every axis. vth = thermal_velocity() = sqrt((Pxx + Pyy + Pzz) / (3 n)), P_dd the second moment
 *   about u along axis d (from the 1-D projections).
 *
 * CANONICAL COEFFICIENT ORDER (on-disk convention, indices are not stored)
 *   Total-degree triangle l + m + n < order_used. Shells s = 0 .. order_used-1 ascending;
 *   inside a shell: l = 0 .. s ascending, then m = 0 .. s-l ascending, n = s-l-m.
 *   The coefficient c[l,m,n] multiplies Hz_l(z) Hy_m(y) Hx_n(x), so coeffs.size() == n_coeffs(order_used).
 *
 * STOPPING RULE
 *   After each shell delta = sqrt(max(0, 1 - sum c^2 / ||f||^2)) (Parseval); the transform exits at the
 *   first shell s > 2 with delta < tolerance (order_used = s + 1, so order_used >= 4), otherwise
 *   order_used = max_order. deltas holds the delta after every computed shell.
 *
 * FAILURE BEHAVIOUR
 *   Nothing throws. adaptive_transform() returns an empty Spectrum (order_used 0, no coeffs, no deltas)
 *   for max_order <= 0, an invalid grid (non-positive shape), a cube whose size differs from the grid,
 *   vth <= 0 or NaN, or a cube with zero power. drift_velocity()/thermal_velocity() return 0 for an
 *   invalid grid or a cube without mass. reconstruct() returns a zero cube for an empty or inconsistent
 *   spectrum (coeffs.size() != n_coeffs(order_used)) or vth <= 0, and an empty vector for an invalid grid.
 *   cube_from_blocks() and global_cell_bbox() return false and leave their outputs untouched.
 *
 * GLOBAL CELL RANGE (file format version 2, columns 5..10 of HERMITEMETA)
 *   global_cell_bbox() maps a Grid to the cell range [lo, hi) on the velocity mesh (centre of global cell i =
 *   mesh_min + (i + 0.5) * mesh_dv). For the grid of cube_from_blocks() this is the block-aligned bounding box of the
 *   present blocks; the writer stores it (as vxmin, vxmax, vymin, vymax, vzmin, vzmax, see bbox_to_meta_range()) so a
 *   reader can zero the reconstruction outside it.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace ASTERIX {
namespace hermite {

struct Spectrum {
   int order_used = 0;          // number of shells kept
   std::vector<double> coeffs;  // canonical order, size n_coeffs(order_used)
   std::vector<double> deltas;  // Parseval delta after each computed shell
};

// Arbitrary Cartesian grid of cell centres (shape may differ per axis).
struct Grid {
   std::array<double, 3> v0;  // centre of the first cell along x, y, z [m/s]
   std::array<double, 3> dv;  // cell size along x, y, z [m/s]
   std::array<int, 3> shape;  // number of cells nx, ny, nz
};

// Number of coefficients in shells 0..order-1: order*(order+1)*(order+2)/6 (0 for order <= 0).
std::size_t n_coeffs(int order);

// ---- Grid-based API ---------------------------------------------------------

// Bulk velocity [ux,uy,uz] (zero vector if the cube has no mass).
std::array<double, 3> drift_velocity(const std::vector<double>& cube, const Grid& grid);

// Isotropic thermal velocity sqrt((Pxx+Pyy+Pzz)/(3n)) about the drift u (0 if the cube has no mass).
double thermal_velocity(const std::vector<double>& cube, const Grid& grid, const std::array<double, 3>& u);

// Level-by-level transform with Parseval stopping (see above).
Spectrum adaptive_transform(const std::vector<double>& cube, const Grid& grid, double vth,
                            const std::array<double, 3>& u, int max_order, double tolerance);

// Dense cube on `grid` from a spectrum. Empty spectrum -> zero cube.
std::vector<double> reconstruct(const Spectrum& spectrum, const Grid& grid, double vth,
                                const std::array<double, 3>& u);

// ---- Cubic API: [-vlim, vlim]^3 sampled with vlen^3 cells (thin wrappers over the Grid API) ----

// Cell-centre velocity axis, size vlen.
std::vector<double> velocity_axis(double vlim, int vlen);

std::array<double, 3> drift_velocity(const std::vector<double>& cube, double vlim, int vlen);
double thermal_velocity(const std::vector<double>& cube, double vlim, int vlen, const std::array<double, 3>& u);
Spectrum adaptive_transform(const std::vector<double>& cube, double vlim, int vlen, double vth,
                            const std::array<double, 3>& u, int max_order, double tolerance);
std::vector<double> reconstruct(const Spectrum& spectrum, double vlim, int vlen, double vth,
                                const std::array<double, 3>& u);

// Assemble a dense cube from velocity-space blocks as stored by Vlasiator.
//   block_params: 6 values per block {vx, vy, vz of the block's min corner, dvx, dvy, dvz}
//                 (dv is the CELL size, as in BlockParams::DVX..DVZ);
//   values:       wid^3 values per block, index k*wid^2 + j*wid + i (x fastest);
//   n_blocks:     number of blocks.
// The grid is the tight bounding box of the present blocks (no margin; absent
// cells are 0). Block order is irrelevant. Returns false (outputs untouched) if
// there are no blocks, wid < 1, the cell size is not identical in all blocks
// (relative 1e-6), or a block corner is not on the lattice of the first block.
template <class Real, class Value>
bool cube_from_blocks(int wid, std::size_t n_blocks, const Real* block_params, const Value* values,
                      std::vector<double>& cube, Grid& grid) {
   if (n_blocks == 0 || wid < 1 || block_params == nullptr || values == nullptr) return false;
   std::array<double, 3> dv, vmin;
   for (int d = 0; d < 3; ++d) {
      dv[d] = static_cast<double>(block_params[3 + d]);
      vmin[d] = static_cast<double>(block_params[d]);
      if (!(dv[d] > 0.0)) return false;
   }
   for (std::size_t b = 1; b < n_blocks; ++b)
      for (int d = 0; d < 3; ++d) {
         if (std::fabs(static_cast<double>(block_params[b * 6 + 3 + d]) - dv[d]) > 1e-6 * dv[d]) return false;
         vmin[d] = std::min(vmin[d], static_cast<double>(block_params[b * 6 + d]));
      }
   // Block origin in cells relative to the lowest corner; reject off-lattice blocks.
   std::vector<std::array<std::size_t, 3>> org(n_blocks);
   std::array<std::size_t, 3> n{0, 0, 0};
   for (std::size_t b = 0; b < n_blocks; ++b)
      for (int d = 0; d < 3; ++d) {
         const double c = (static_cast<double>(block_params[b * 6 + d]) - vmin[d]) / dv[d];
         const double r = std::round(c);
         if (std::fabs(c - r) > 1e-3) return false;
         org[b][d] = static_cast<std::size_t>(r);
         n[d] = std::max(n[d], org[b][d] + static_cast<std::size_t>(wid));
      }
   cube.assign(n[0] * n[1] * n[2], 0.0);
   const std::size_t W = wid;
   for (std::size_t b = 0; b < n_blocks; ++b)
      for (std::size_t k = 0; k < W; ++k)
         for (std::size_t j = 0; j < W; ++j)
            for (std::size_t i = 0; i < W; ++i)
               cube[((org[b][2] + k) * n[1] + org[b][1] + j) * n[0] + org[b][0] + i] =
                  static_cast<double>(values[b * W * W * W + k * W * W + j * W + i]);
   for (int d = 0; d < 3; ++d) {
      grid.v0[d] = vmin[d] + 0.5 * dv[d];
      grid.dv[d] = dv[d];
      grid.shape[d] = static_cast<int>(n[d]);
   }
   return true;
}

// Global velocity-mesh cell range covered by `grid`.
//   mesh_min: lower corner of the velocity mesh (meshMinLimits); mesh_dv: mesh cell size.
//   bbox = {lo_x, lo_y, lo_z, hi_x, hi_y, hi_z}: lower bound inclusive, upper bound EXCLUSIVE, in global cell indices.
// Index convention: global cell i along an axis has its centre at mesh_min + (i + 0.5) * mesh_dv, so
//   lo = (grid.v0 - 0.5 * dv - mesh_min) / mesh_dv and hi = lo + shape. When the grid comes from cube_from_blocks()
// the range is block aligned (multiples of WID) because it is the bounding box of whole blocks.
// Returns false (bbox untouched) for an empty grid, a cell size that differs from mesh_dv (relative 1e-6), a grid off
// the mesh lattice (tolerance 1e-3 cell) or one that starts below the mesh origin.
inline bool global_cell_bbox(const Grid& grid, const std::array<double, 3>& mesh_min,
                             const std::array<double, 3>& mesh_dv, std::array<std::uint32_t, 6>& bbox) {
   std::array<std::uint32_t, 6> out{};
   for (int d = 0; d < 3; ++d) {
      if (grid.shape[d] <= 0 || !(mesh_dv[d] > 0.0)) return false;
      if (!(std::fabs(grid.dv[d] - mesh_dv[d]) <= 1e-6 * mesh_dv[d])) return false;
      const double c = (grid.v0[d] - 0.5 * mesh_dv[d] - mesh_min[d]) / mesh_dv[d];
      const double r = std::round(c);
      if (!(std::fabs(c - r) <= 1e-3) || r < 0.0 || r > 4.0e9) return false;
      out[d] = static_cast<std::uint32_t>(r);
      out[d + 3] = out[d] + static_cast<std::uint32_t>(grid.shape[d]);
   }
   bbox = out;
   return true;
}

// Reorders a global_cell_bbox() result {lo_x, lo_y, lo_z, hi_x, hi_y, hi_z} into the HERMITEMETA column order
// {vxmin, vxmax, vymin, vymax, vzmin, vzmax} (min inclusive, max EXCLUSIVE, global velocity-mesh cell indices) as float32.
// This is the single place where that reordering happens. float32 represents integers exactly up to 2^24; returns
// false (range untouched) if any index is larger.
inline bool bbox_to_meta_range(const std::array<std::uint32_t, 6>& bbox, std::array<float, 6>& range) {
   constexpr std::uint32_t max_exact = 1u << 24;
   for (std::uint32_t v : bbox)
      if (v > max_exact) return false;
   for (int d = 0; d < 3; ++d) {
      range[2 * d] = static_cast<float>(bbox[d]);         // vXmin = lo
      range[2 * d + 1] = static_cast<float>(bbox[d + 3]); // vXmax = hi (exclusive)
   }
   return true;
}

} // namespace hermite
} // namespace ASTERIX
