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
 *   c = dv_x dv_y dv_z * sum g * phi_l(z) phi_m(y) phi_n(x), g = sqrt(f) (see EXPANSION), i.e. in units of
 *   sqrt(f) * (m/s)^(3/2). Hence a cube that is zero outside a sub-box may be cropped to that box without changing any
 *   result.
 *   The cubic API taking (vlim, vlen) is a thin wrapper over
 *   Grid{v0 = -vlim + dv/2, dv = 2*vlim/vlen, shape = vlen^3}.
 *
 * BASIS
 *   One isotropic thermal velocity vth (the standard deviation of the Gaussian) and the bulk velocity u
 *   define the 1-D orthonormal Hermite functions phi_n((v - u_d) / s) / sqrt(s), n = 0, 1, ..., used along every
 *   axis, with the scale s = sqrt(2) * vth (the width of the square root of a Maxwellian). vth = thermal_velocity() =
 *   sqrt((Pxx + Pyy + Pzz) / (3 n)), P_dd the second moment about u along axis d (from the 1-D projections); the
 *   functions take the physical vth of f and apply the factor sqrt(2) themselves.
 *
 * EXPANSION
 *   The series is that of g = sqrt(max(f, 0)) (port of the sqrt(f) transform of faiser vdf_tools.py); reconstruct() returns
 *   g^2. The reconstruction is non-negative by construction, its ringing outside the support is ~eps^2 instead of ~eps
 *   (eps: the relative size of the ringing of an expansion of f itself), and delta (Parseval, ||g||^2 = n) is the L2
 *   error of the g series, i.e. a density misfit; the L2 error in f is larger (1.6-2.4x in faiser tests).
 *
 * CANONICAL COEFFICIENT ORDER (on-disk convention, indices are not stored)
 *   Total-degree triangle l + m + n < order_used. Shells s = 0 .. order_used-1 ascending;
 *   inside a shell: l = 0 .. s ascending, then m = 0 .. s-l ascending, n = s-l-m.
 *   The coefficient c[l,m,n] multiplies Hz_l(z) Hy_m(y) Hx_n(x), so coeffs.size() == n_coeffs(order_used).
 *
 * STOPPING RULE
 *   After each shell delta = sqrt(max(0, 1 - sum c^2 / ||g||^2)) (Parseval); the transform exits at the
 *   first shell s > 2 with delta < tolerance (order_used = s + 1, so order_used >= 4), otherwise
 *   order_used = max_order. deltas holds the delta after every computed shell.
 *
 * FAILURE BEHAVIOUR
 *   Nothing throws. adaptive_transform() returns an empty Spectrum (order_used 0, no coeffs, no deltas)
 *   for max_order <= 0, an invalid grid (non-positive shape), a cube whose size differs from the grid,
 *   vth <= 0 or NaN, or a cube with zero power. drift_velocity()/thermal_velocity() return 0 for an
 *   invalid grid or a cube without mass. reconstruct() returns a zero cube for an empty or inconsistent
 *   spectrum (coeffs.size() != n_coeffs(order_used)) or vth <= 0, and an empty vector for an invalid grid.
 *   cube_from_blocks() returns false and leaves its outputs untouched.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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

// Level-by-level transform of g = sqrt(max(f, 0)) with Parseval stopping (see above); `vth` is the physical thermal
// velocity of f.
Spectrum adaptive_transform(const std::vector<double>& cube, const Grid& grid, double vth,
                            const std::array<double, 3>& u, int max_order, double tolerance);

// Dense cube of f = g^2 >= 0 on `grid` from a spectrum. Empty spectrum -> zero cube.
std::vector<double> reconstruct(const Spectrum& spectrum, const Grid& grid, double vth,
                                const std::array<double, 3>& u);

// Cesaro (C,1) = Fejer mean of the series: coefficient of shell s = l + m + n is weighted by (N - s) / N, N = order.
// `coeffs` is in the canonical order with n_coeffs(order) entries (size mismatch -> returned unchanged).
std::vector<double> cesaro_coeffs(const std::vector<double>& coeffs, int order);

// reconstruct() for several coefficient sets that share order, vth, u and grid with `base` (its own coeffs are not
// used): the basis functions are computed once. Result k belongs to coeff_sets[k]; a set of the wrong size gives a
// zero cube.
std::vector<std::vector<double>> reconstruct_sets(const Spectrum& base, const std::vector<std::vector<double>>& coeff_sets,
                                                  const Grid& grid, double vth, const std::array<double, 3>& u);

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
// The cube is cropped to the present blocks (no margin; absent cells inside the
// crop are 0), grid describes its cell centres. Block order is irrelevant. Returns false (outputs untouched) if
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

} // namespace hermite
} // namespace ASTERIX
