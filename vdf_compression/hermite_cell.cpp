/*
 * Spatial-cell glue for the adaptive Hermite transform; see hermite_cell.h.
 *
 * The VDF cube is built directly from the velocity blocks (hermite::cube_from_blocks) rather than with
 * extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(): the cube is cropped to the blocks that are present
 * (the coefficients are the same as for the whole mesh, the VDF being zero elsewhere, but it is much cheaper), and
 * the extractor reads dv from block index 1 (compression_tools.cpp), which is out of bounds for single-block cells.
 *
 * Geometry: BlockParams::VXCRD..VZCRD is the min corner of a block and BlockParams::DVX..DVZ is the CELL size,
 * so cell (i,j,k) of a block has its centre at corner + (i + 0.5) * dv.
 */
#include "hermite_cell.h"

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <iostream>

#include "../common.h"
#include "../object_wrapper.h"
#include "../spatial_cells/spatial_cell_wrapper.hpp"
#include "../spatial_cells/velocity_block_container.h"
#include "../velocity_mesh_parameters.h"
#include "hermite.h"

namespace ASTERIX {

std::optional<HermiteCell> compute_hermite_cell(spatial_cell::SpatialCell* sc, unsigned int popID, int max_order,
                                                double tolerance) {
   if (sc == nullptr || max_order <= 0) return std::nullopt;
   if (sc->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) return std::nullopt;

   const vmesh::VelocityBlockContainer* blocks = sc->get_velocity_blocks(popID);
   const std::size_t n_blocks = blocks->size();
   if (n_blocks == 0) return std::nullopt;

   static_assert(BlockParams::N_VELOCITY_BLOCK_PARAMS == 6, "cube_from_blocks expects 6 parameters per block");
   std::vector<double> cube;
   hermite::Grid grid;
   if (!hermite::cube_from_blocks<Real, Realf>(WID, n_blocks, sc->get_block_parameters(popID), blocks->getData(),
                                               cube, grid))
      return std::nullopt;

   const auto u = hermite::drift_velocity(cube, grid);
   const double vth = hermite::thermal_velocity(cube, grid, u);
   if (!(vth > 0.0)) return std::nullopt;

   const hermite::Spectrum s = hermite::adaptive_transform(cube, grid, vth, u, max_order, tolerance);
   if (s.order_used <= 0 || s.coeffs.empty()) return std::nullopt;

   HermiteCell out;
   for (int d = 0; d < 3; ++d) out.meta[d] = static_cast<float>(u[d]);
   out.meta[3] = static_cast<float>(vth);
   out.meta[4] = static_cast<float>(s.deltas.back());
   out.order_used = s.order_used;
   out.coeffs.assign(s.coeffs.begin(), s.coeffs.end());
   assert(out.coeffs.size() == hermite::n_coeffs(out.order_used));
   return out;
}

// ---------------------------------------------------------------------------------------------------------------------
// Round trip
// ---------------------------------------------------------------------------------------------------------------------

#ifndef USE_GPU
namespace {

// The velocity mesh as a grid of cell centres. A cube on it is indexed cube[(iz * ny + iy) * nx + ix] (= global cell
// index) and is cut into WID^3 blocks, block (bx,by,bz) being number (bz * nby + by) * nbx + bx.
using MeshGrid = hermite::Grid;

struct BlockCube {
   std::size_t nx, ny, nbx, nby, nbz;
   explicit BlockCube(const MeshGrid& g)
      : nx(g.shape[0]), ny(g.shape[1]), nbx(g.shape[0] / WID), nby(g.shape[1] / WID), nbz(g.shape[2] / WID) {}
   std::size_t n_blocks() const { return nbx * nby * nbz; }
   std::size_t block(std::size_t bx, std::size_t by, std::size_t bz) const { return (bz * nby + by) * nbx + bx; }
   // Index in the cube of cell (i,j,k) of block (bx,by,bz).
   std::size_t cell(std::size_t bx, std::size_t by, std::size_t bz, std::size_t i, std::size_t j, std::size_t k) const {
      return ((bz * WID + k) * ny + by * WID + j) * nx + bx * WID + i;
   }
   template <class F> void for_each_block(F f) const {
      for (std::size_t bz = 0; bz < nbz; ++bz)
         for (std::size_t by = 0; by < nby; ++by)
            for (std::size_t bx = 0; bx < nbx; ++bx) f(bx, by, bz);
   }
};

// Blocks with at least one cell on or above the threshold.
std::vector<char> blocks_with_content(const BlockCube& b, const std::vector<double>& f, double threshold) {
   std::vector<char> content(b.n_blocks(), 0);
   b.for_each_block([&](std::size_t bx, std::size_t by, std::size_t bz) {
      for (std::size_t k = 0; k < WID; ++k)
         for (std::size_t j = 0; j < WID; ++j)
            for (std::size_t i = 0; i < WID; ++i)
               if (f[b.cell(bx, by, bz, i, j, k)] >= threshold) {
                  content[b.block(bx, by, bz)] = 1;
                  return;
               }
   });
   return content;
}

// All blocks within `width` blocks (in all 26 directions) of a block of `content`.
std::vector<char> dilate(const BlockCube& b, const std::vector<char>& content, int width) {
   std::vector<char> out(content.size(), 0);
   const long w = std::max(width, 0);
   b.for_each_block([&](std::size_t bx, std::size_t by, std::size_t bz) {
      if (!content[b.block(bx, by, bz)]) return;
      for (long dz = -w; dz <= w; ++dz)
         for (long dy = -w; dy <= w; ++dy)
            for (long dx = -w; dx <= w; ++dx) {
               const long cx = static_cast<long>(bx) + dx, cy = static_cast<long>(by) + dy, cz = static_cast<long>(bz) + dz;
               if (cx < 0 || cy < 0 || cz < 0 || cx >= static_cast<long>(b.nbx) || cy >= static_cast<long>(b.nby) ||
                   cz >= static_cast<long>(b.nbz))
                  continue;
               out[b.block(cx, cy, cz)] = 1;
            }
   });
   return out;
}

// The blocks of the hybrid reconstruction: their indices on the velocity mesh and their data.
struct KeptBlocks {
   std::vector<std::array<vmesh::GlobalID, 3>> index;
   std::vector<Realf> data; // WID^3 values per block, cell (i,j,k) at k*WID^2 + j*WID + i
};

// Hybrid reconstruction (see hermite_roundtrip in hermite_cell.h): the Cesaro mean gives the sparsity contour and the
// block mask; core (inside the contour) = plain series, halo (inside the mask, outside the contour) = Cesaro values.
bool hybrid_blocks(const hermite::Spectrum& spectrum, const MeshGrid& mesh, const std::array<double, 3>& u, double vth,
                   double threshold, int width, KeptBlocks& out) {
   if (mesh.shape[0] % WID || mesh.shape[1] % WID || mesh.shape[2] % WID) return false; // the mesh is made of blocks
   const BlockCube b(mesh);

   // Plain and Cesaro series share the basis functions.
   const auto sets = hermite::reconstruct_sets(
      spectrum, {spectrum.coeffs, hermite::cesaro_coeffs(spectrum.coeffs, spectrum.order_used)}, mesh, vth, u);
   const std::size_t n_cells = b.nx * b.ny * b.nbz * WID;
   if (sets.size() != 2 || sets[0].size() != n_cells || sets[1].size() != n_cells) return false;
   const std::vector<double>& rec = sets[0];
   const std::vector<double>& ces = sets[1];

   const std::vector<char> keep = dilate(b, blocks_with_content(b, ces, threshold), width);
   b.for_each_block([&](std::size_t bx, std::size_t by, std::size_t bz) {
      if (!keep[b.block(bx, by, bz)]) return;
      out.index.push_back({static_cast<vmesh::GlobalID>(bx), static_cast<vmesh::GlobalID>(by), static_cast<vmesh::GlobalID>(bz)});
      for (std::size_t k = 0; k < WID; ++k)
         for (std::size_t j = 0; j < WID; ++j)
            for (std::size_t i = 0; i < WID; ++i) {
               const std::size_t c = b.cell(bx, by, bz, i, j, k);
               out.data.push_back(static_cast<Realf>(ces[c] >= threshold ? rec[c] : ces[c]));
            }
   });
   return true;
}

// Round trip of one simulated cell. Returns nullptr on success, otherwise the reason why the cell could not be compressed.
const char* roundtrip_cell(spatial_cell::SpatialCell* sc, unsigned int popID, int max_order, double tolerance,
                           double threshold, int width, const MeshGrid& mesh, RoundTripStats& st) {
   const std::size_t blocksBefore = sc->get_number_of_velocity_blocks(popID);
   if (blocksBefore == 0) return "the cell has no velocity blocks";

   // Forward transform exactly as for the file: what is reconstructed is what a file would hold (float32 coefficients).
   const std::optional<HermiteCell> cell = compute_hermite_cell(sc, popID, max_order, tolerance);
   if (!cell) return "the transform failed (no thermal spread, or velocity blocks that are not on a uniform lattice)";
   hermite::Spectrum spectrum;
   spectrum.order_used = cell->order_used;
   spectrum.coeffs.assign(cell->coeffs.begin(), cell->coeffs.end());
   const std::array<double, 3> u{cell->meta[0], cell->meta[1], cell->meta[2]};

   KeptBlocks kept;
   if (!hybrid_blocks(spectrum, mesh, u, cell->meta[3], threshold, width, kept)) return "the hybrid reconstruction failed";
   std::vector<vmesh::GlobalID> gids;
   for (auto& idx : kept.index) gids.push_back(sc->get_velocity_block(popID, idx.data()));
   if (gids.empty()) return "no velocity block is left after the reconstruction";

   // Replace the blocks of the cell by the block mask. The caller must run adjustVelocityBlocks() (grid.h) afterwards: it
   // applies the sparsity rule with the spatial neighbours and re-synchronizes the block lists of the ghost copies on
   // the other ranks, which the next translation relies on.
   sc->clear(popID);
   sc->add_velocity_blocks<Realf>(popID, gids, kept.data.data());

   ++st.cells;
   st.blocks_before += blocksBefore;
   st.sum_order += cell->order_used;
   return nullptr;
}

} // namespace
#endif

RoundTripStats hermite_roundtrip(const std::vector<std::uint64_t>& cellIds,
                                 const std::vector<spatial_cell::SpatialCell*>& cells, unsigned int popID, int max_order,
                                 double tolerance) {
   RoundTripStats total;
#ifdef USE_GPU
   // Parameter parsing aborts before this point in GPU builds; kept as a safe no-op.
   (void)cellIds;
   (void)cells;
   (void)popID;
   (void)max_order;
   (void)tolerance;
   return total;
#else
   const auto& species = getObjectWrapper().particleSpecies[popID];
   const double threshold = static_cast<double>(species.sparseMinValue);
   const int width = species.sparseBlockAddWidthV;
   const auto& vmeshInfo = vmesh::getMeshWrapper()->velocityMeshes->at(species.velocityMesh);
   MeshGrid mesh;
   for (int d = 0; d < 3; ++d) {
      mesh.v0[d] = vmeshInfo.meshMinLimits[d] + 0.5 * vmeshInfo.cellSize[d];
      mesh.dv[d] = vmeshInfo.cellSize[d];
      mesh.shape[d] = static_cast<int>(vmeshInfo.gridLength[d] * vmeshInfo.blockLength[d]);
   }

   // Every simulated cell is round-tripped; the cells are independent (each one only touches its own blocks).
   // DO_NOT_COMPUTE cells have no velocity space and are not part of the simulation.
   std::size_t nFailed = 0;
   std::uint64_t firstFailedCell = 0;
   const char* firstReason = nullptr;
   #pragma omp parallel
   {
      RoundTripStats local;
      #pragma omp for schedule(dynamic)
      for (std::size_t c = 0; c < cells.size(); ++c) {
         if (cells[c]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) continue;
         const char* reason = roundtrip_cell(cells[c], popID, max_order, tolerance, threshold, width, mesh, local);
         if (reason) {
            #pragma omp critical(hermite_roundtrip_failure)
            if (nFailed++ == 0) {
               firstFailedCell = cellIds[c];
               firstReason = reason;
            }
         }
      }
      #pragma omp critical(hermite_roundtrip_stats)
      {
         total.cells += local.cells;
         total.blocks_before += local.blocks_before;
         total.sum_order += local.sum_order;
      }
   }

   // A simulated cell that cannot be compressed is fatal: the run would no longer be a Hermite-compressed run.
   if (nFailed > 0) {
      std::cerr << "ERROR: Hermite round trip cannot compress " << nFailed << " cell(s) on this rank, the first is cell "
                << firstFailedCell << ": " << firstReason << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   // Sum over the MPI ranks (all values are exact in a double).
   double sums[3] = {static_cast<double>(total.cells), static_cast<double>(total.blocks_before), total.sum_order};
   double all[3];
   MPI_Allreduce(sums, all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   total.cells = static_cast<std::uint64_t>(all[0]);
   total.blocks_before = static_cast<std::uint64_t>(all[1]);
   total.sum_order = all[2];
   return total;
#endif
}

} // namespace ASTERIX
