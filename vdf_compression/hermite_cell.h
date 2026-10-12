#pragma once
/*
 * Glue between the pure Hermite transform (hermite.h) and a Vlasiator
 * SpatialCell. This is the only Hermite file that depends on Vlasiator types.
 *
 * File format history (HERMITECOEFFS attribute "version"):
 *   1: HERMITEMETA = float32 x5  (ux, uy, uz, vth, delta)
 *   2: HERMITEMETA = float32 x11 (ux, uy, uz, vth, delta, vxmin, vxmax, vymin, vymax, vzmin, vzmax); the last six are
 *      the velocity-space bounding box in global velocity-mesh cell indices (min inclusive, max exclusive).
 *      (A short-lived intermediate layout stored them in a separate HERMITEBBOX array; it is no longer produced.)
 */

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace spatial_cell {
class SpatialCell;
}

namespace ASTERIX {

struct HermiteCell {
   // Same order as the 11 columns of HERMITEMETA (format version 2):
   //   ux, uy, uz       [m/s] bulk velocity
   //   vth              [m/s] isotropic thermal velocity sqrt((Pxx+Pyy+Pzz)/(3n))
   //   delta            Parseval delta sqrt(max(0, 1 - sum c^2 / ||f||^2)) at order_used
   //   vxmin, vxmax, vymin, vymax, vzmin, vzmax
   //                    bounding box of the velocity blocks present in the cell, in GLOBAL velocity-mesh CELL indices
   //                    (stored as float32, exact up to 2^24); *min inclusive, *max EXCLUSIVE; block aligned (multiples
   //                    of WID). Global cell i has its centre at meshMinLimits + (i + 0.5) * cellSize. Equals the cell
   //                    range of the cube the transform ran on.
   std::array<float, 11> meta{0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
   int order_used = 0;        // number of shells kept
   std::vector<float> coeffs; // canonical order (see hermite.h), size n_coeffs(order_used)
};

// Adaptive Hermite transform of population `popID` of `sc`.
// Returns an empty optional for DO_NOT_COMPUTE cells, cells without velocity blocks,
// cells with zero power, non-uniform / non-lattice velocity blocks or blocks off the velocity-mesh lattice.
std::optional<HermiteCell> compute_hermite_cell(spatial_cell::SpatialCell* sc, unsigned int popID, int max_order,
                                                double tolerance);

} // namespace ASTERIX
