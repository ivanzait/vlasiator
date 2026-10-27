#pragma once
/*
 * Glue between the pure Hermite transform (hermite.h) and a Vlasiator
 * SpatialCell. This is the only Hermite file that depends on Vlasiator types.
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
   // The 5 columns of HERMITEMETA:
   //   ux, uy, uz       [m/s] bulk velocity
   //   vth              [m/s] isotropic thermal velocity sqrt((Pxx+Pyy+Pzz)/(3n)) of f
   //   delta            Parseval delta sqrt(max(0, 1 - sum c^2 / ||g||^2)) at order_used, g = sqrt(f)
   std::array<float, 5> meta{0.f, 0.f, 0.f, 0.f, 1.f};
   int order_used = 0;        // number of shells kept
   std::vector<float> coeffs; // canonical order (see hermite.h), size n_coeffs(order_used)
};

// Adaptive Hermite transform of population `popID` of `sc`.
// Returns an empty optional for DO_NOT_COMPUTE cells, cells without velocity blocks, cells with zero power or
// non-uniform / non-lattice velocity blocks.
std::optional<HermiteCell> compute_hermite_cell(spatial_cell::SpatialCell* sc, unsigned int popID, int max_order,
                                                double tolerance);

// Totals of one round-trip pass (summed over all MPI ranks).
struct RoundTripStats {
   std::uint64_t cells = 0;         // cells that were round-tripped (every simulated cell)
   std::uint64_t blocks_before = 0; // velocity blocks of those cells before the round trip
   double sum_order = 0.0;          // sum of the Hermite orders used
};

// Hermite round trip of population `popID` in every simulated cell of `cells` (cellIds are the corresponding cell ids, only
// used in error messages; DO_NOT_COMPUTE cells have no velocity space and are left out). OpenMP over cells, collective MPI
// reduction of the statistics over MPI_COMM_WORLD: all ranks must call it. A simulated cell that cannot be compressed is a
// fatal error (message with the cell id and the reason, then MPI_Abort). Per cell: forward adaptive transform (max_order,
// tolerance), then the hybrid reconstruction:
//   - the Cesaro (C,1) mean of the series defines the sparsity contour (Cesaro value >= the population's
//     sparseMinValue) and the block mask: a block has content if any of its cells is inside the contour, and every
//     block within the population's sparse blockAddWidthV blocks of a content block is kept whole (all 26 directions);
//   - inside the mask the core (inside the contour) is reconstructed with the plain series and the halo (inside the
//     mask, outside the contour) with the Cesaro values; everything outside the mask is zero.
// Both series are evaluated on the whole velocity mesh; the mask sets the boundary. The cell's blocks are replaced by the
// blocks of the mask. The forward transform is compute_hermite_cell(), so what is reconstructed is what a file holds. The caller must then run adjustVelocityBlocks(mpiGrid, cells, true, popID) (grid.h), which applies
// Vlasiator's sparsity rule including the spatial neighbours and re-synchronizes the block lists of the ghost copies on
// the other ranks (without it the next translation can abort in the block exchange). Moments are not rescaled. Not
// supported in GPU builds.
RoundTripStats hermite_roundtrip(const std::vector<std::uint64_t>& cellIds,
                                 const std::vector<spatial_cell::SpatialCell*>& cells, unsigned int popID, int max_order,
                                 double tolerance);

} // namespace ASTERIX
