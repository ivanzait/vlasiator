/*
 * Spatial-cell glue for the adaptive Hermite transform; see hermite_cell.h.
 *
 * The VDF cube is built directly from the velocity blocks (hermite::cube_from_blocks) rather than with
 * extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(): this gives exact integer cell indices on the
 * block lattice and the tight bounding box of the present blocks, without margin. The extractor also reads
 * dv from block index 1 (compression_tools.cpp), which is out of bounds for single-block cells.
 *
 * Geometry: BlockParams::VXCRD..VZCRD is the min corner of a block and BlockParams::DVX..DVZ is the CELL size,
 * so cell (i,j,k) of a block has its centre at corner + (i + 0.5) * dv.
 *
 * The cube's global cell range on the velocity mesh is stored in HermiteCell::meta[5..10] (format version 2).
 */
#include "hermite_cell.h"

#include <cassert>

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

   // Global cell range of the cube on the velocity mesh (block aligned; see hermite::global_cell_bbox).
   const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
   const auto& mesh = vmesh::getMeshWrapper()->velocityMeshes->at(meshID);
   const std::array<double, 3> meshMin{mesh.meshMinLimits[0], mesh.meshMinLimits[1], mesh.meshMinLimits[2]};
   const std::array<double, 3> meshDv{mesh.cellSize[0], mesh.cellSize[1], mesh.cellSize[2]};
   std::array<uint32_t, 6> bbox;
   if (!hermite::global_cell_bbox(grid, meshMin, meshDv, bbox)) return std::nullopt;
   std::array<float, 6> range; // vxmin, vxmax, vymin, vymax, vzmin, vzmax
   if (!hermite::bbox_to_meta_range(bbox, range)) return std::nullopt; // indices do not fit float32 exactly

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
   for (int i = 0; i < 6; ++i) out.meta[5 + i] = range[i];
   out.coeffs.assign(s.coeffs.begin(), s.coeffs.end());
   assert(out.coeffs.size() == hermite::n_coeffs(out.order_used));
   return out;
}

} // namespace ASTERIX
