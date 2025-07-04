#include "vdf_tools.h"
#include "../object_wrapper.h"
#include <array>
#include <cmath> 


////////////// Vector Algebra
using Vector3 = std::array<float, 3>;
using Matrix3x3 = std::array<Vector3, 3>; 

// linspace
std::vector<float> linspace(float start, float end, int len) {
   const float step = (end - start) / (len - 1);
   std::vector<float> x(len,0.0f);
   for (std::size_t i = 0; i < x.size(); ++i) {
      x[i] = start + i * step;
   }
   return x;
}

Vector3 normalize(const Vector3& v) {
    float norm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return {v[0]/norm, v[1]/norm, v[2]/norm};
}

Vector3 cross(const Vector3& a, const Vector3& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

Matrix3x3 build_rotation_matrix(const Vector3& E, const Vector3& B) {
    Vector3 b_hat = normalize(B);
    Vector3 e_hat = normalize(E);
    Vector3 exb = cross(e_hat, b_hat);
    Vector3 exb_hat = normalize(exb);
    Vector3 perp = cross(b_hat, exb_hat); 

    return { b_hat, exb_hat, perp }; 
}


/// Matrix transpose
Matrix3x3 transpose(const Matrix3x3& R) {
    Matrix3x3 Rt;
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
        Rt[i][j] = R[j][i];
    return Rt;
}

// Helper to get linear index from 3D index
size_t idx(const OrderedVDF& vdf, size_t i, size_t j, size_t k) {
    return i * (vdf.shape[1] * vdf.shape[2]) + j * vdf.shape[2] + k;
}

// Trilinear interpolation in original VDF at point (vx, vy, vz)
float trilinear_interpolate(const OrderedVDF& vdf, float vx, float vy, float vz) {
    // Normalize vx,vy,vz to grid indices:
    float x_frac = (vx - vdf.v_limits[0]) / (vdf.v_limits[3] - vdf.v_limits[0]) * (vdf.shape[0] - 1);
    float y_frac = (vy - vdf.v_limits[1]) / (vdf.v_limits[4] - vdf.v_limits[1]) * (vdf.shape[1] - 1);
    float z_frac = (vz - vdf.v_limits[2]) / (vdf.v_limits[5] - vdf.v_limits[2]) * (vdf.shape[2] - 1);

    int i0 = std::floor(x_frac);
    int j0 = std::floor(y_frac);
    int k0 = std::floor(z_frac);

    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;

    // Check bounds
    if (i0 < 0 || j0 < 0 || k0 < 0 || i1 >= (int)vdf.shape[0] || j1 >= (int)vdf.shape[1] || k1 >= (int)vdf.shape[2])
        return 0.0f;  // Outside grid → return 0 or some fallback

    float xd = x_frac - i0;
    float yd = y_frac - j0;
    float zd = z_frac - k0;

    // Fetch the 8 corners
    float c000 = vdf.vdf_vals[idx(vdf, i0, j0, k0)];
    float c100 = vdf.vdf_vals[idx(vdf, i1, j0, k0)];
    float c010 = vdf.vdf_vals[idx(vdf, i0, j1, k0)];
    float c110 = vdf.vdf_vals[idx(vdf, i1, j1, k0)];
    float c001 = vdf.vdf_vals[idx(vdf, i0, j0, k1)];
    float c101 = vdf.vdf_vals[idx(vdf, i1, j0, k1)];
    float c011 = vdf.vdf_vals[idx(vdf, i0, j1, k1)];
    float c111 = vdf.vdf_vals[idx(vdf, i1, j1, k1)];

    // Interpolate along x
    float c00 = c000 * (1 - xd) + c100 * xd;
    float c10 = c010 * (1 - xd) + c110 * xd;
    float c01 = c001 * (1 - xd) + c101 * xd;
    float c11 = c011 * (1 - xd) + c111 * xd;

    // Interpolate along y
    float c0 = c00 * (1 - yd) + c10 * yd;
    float c1 = c01 * (1 - yd) + c11 * yd;

    // Interpolate along z
    float c = c0 * (1 - zd) + c1 * zd;

    return c;
}



/////////////////// END of vector algebra

//// KOSTIS' EXTRACT
HERMITE::OrderedVDF HERMITE::extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(SpatialCell* sc, uint popID,
                                                                                       int zoom) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   if (zoom != 1) {
      throw std::runtime_error("Zoom is not supported yet!");
   }
   auto* blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer->size();
   const Real* blockParams = sc->get_block_parameters(popID);

   Real sparse = getObjectWrapper().particleSpecies[popID].sparseMinValue;

   // xmin,ymin,zmin,xmax,ymax,zmax;
   std::array<Real, 6> vlims{std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                             std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::lowest(),
                             std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

   // This pass is computing the active vmesh limits
   // Store dvx, dvy, dvz here
   const Real dvx = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVX];
   const Real dvy = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVY];
   const Real dvz = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVZ];
   for (std::size_t n = 0; n < total_blocks; ++n) {
      const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               const Real vx = bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX];
               const Real vy = bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY];
               const Real vz = bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ];
               vlims[0] = std::min(vlims[0], vx);
               vlims[1] = std::min(vlims[1], vy);
               vlims[2] = std::min(vlims[2], vz);
               vlims[3] = std::max(vlims[3], vx);
               vlims[4] = std::max(vlims[4], vy);
               vlims[5] = std::max(vlims[5], vz);
            }
         }
      }
   } // over blocks

   // assert(isPow2(static_cast<size_t>(std::abs(zoom))));
   float ratio = (zoom > 0) ? static_cast<float>(std::abs(zoom)) : 1.0 / static_cast<float>(std::abs(zoom));
   assert(ratio > 0);

   const Real target_dvx = dvx * ratio;
   const Real target_dvy = dvy * ratio;
   const Real target_dvz = dvz * ratio;
   std::size_t nx = std::ceil((vlims[3] - vlims[0]) / target_dvx);
   std::size_t ny = std::ceil((vlims[4] - vlims[1]) / target_dvy);
   std::size_t nz = std::ceil((vlims[5] - vlims[2]) / target_dvz);
   // printf("VDF min box is %zu , %zu %zu \n ", nx, ny, nz);


   Realf* data = blockContainer->getData();
   std::vector<Realf> vspace(nx * ny * nz, Realf(0));
   for (std::size_t n = 0; n < total_blocks; ++n) {
      const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      const Realf* vdf_data = &data[n * WID3];
      const vmesh::GlobalID gid = sc->get_velocity_block_global_id(n, popID);
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               const Real vx = bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX];
               const Real vy = bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY];
               const Real vz = bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ];
               const size_t bbox_i = std::min(static_cast<size_t>(std::floor((vx - vlims[0]) / target_dvx)), nx - 1);
               const size_t bbox_j = std::min(static_cast<size_t>(std::floor((vy - vlims[1]) / target_dvy)), ny - 1);
               const size_t bbox_k = std::min(static_cast<size_t>(std::floor((vz - vlims[2]) / target_dvz)), nz - 1);

               // Averaging
               if (ratio >= 1.0) {
                  const size_t index = bbox_i * (ny * nz) + bbox_j * nz + bbox_k;
                  if (vdf_data[cellIndex(i, j, k)] < sparse) {
                     vspace.at(index) = 0.0;
                  } else {
                     vspace.at(index) += vdf_data[cellIndex(i, j, k)] / ratio;
                  }
                  // vspace.at(index) += vdf_data[cellIndex(i, j, k)] / ratio; // that my sh**
               } else {
                  // Same value in all bins
                  int max_off = 1 / ratio;
                  for (int off_z = 0; off_z <= max_off; off_z++) {
                     for (int off_y = 0; off_y <= max_off; off_y++) {
                        for (int off_x = 0; off_x <= max_off; off_x++) {
                           const size_t index = (bbox_i + off_x) * (ny * nz) + (bbox_j + off_y) * nz + (bbox_k + off_z);
                           if (index < vspace.size()) {
                              if (vdf_data[cellIndex(i, j, k)] < sparse) {
                                 vspace.at(index) = 0.0;
                              } else {
                                 vspace.at(index) =
                                     vdf_data[cellIndex(i, j, k)]; // it was only this thing here without sparcity
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   } // over blocks

   // making log10 !!!
   // for (uint i = 0; i < vspace.size(); ++i) {
   //    vspace[i] = std::log10(std::max(vspace[i], static_cast<float>(sparse))) - std::log10(sparse);
   // }

   return HERMITE::OrderedVDF{.blocks_to_ignore = {},
                              .sparse_vdf_bytes = total_blocks * WID * WID * WID * sizeof(Realf),
                              .vdf_vals = vspace,
                              .v_limits = vlims,
                              .shape = {nx, ny, nz}};
}

void dump_vdf_to_binary_file(const char* filename, uint popID, CellID cid,
                             dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   SpatialCell* sc = mpiGrid[cid];
   assert(sc && "Invalid Pointer to Spatial Cell !");
   HERMITE::OrderedVDF vdf = HERMITE::extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(sc, popID, 1);
   // vdf.save_to_file(filename);
}

////// END of Kostis' extract

////// HERMITE BASIS
std::size_t factorial(std::size_t n) { return (n == 0) ? (1) : (n * factorial(n - 1)); }

// linspace
//std::vector<float> linspace(float start, float end, int len) {
//   const float step = (end - start) / (len - 1);
//   std::vector<float> x(len,0.0f);
//   for (std::size_t i = 0; i < x.size(); ++i) {
//      x[i] = start + i * step;
//   }
//   return x;
//}

// generate hermite polinomials up to the given order
std::vector<std::vector<float>> hermite(std::vector<float>& x, int order) {
   // Recurrence relation: H_n+1 = 2*x*H_n -2*n*H_n-1
   // H0 = 1, H1 = 2x
   // base function with Gauss weights: H_n * exp(-0.5*v^2)
   std::vector<std::vector<float>> hp(order, std::vector<float>(x.size()));
   for (size_t i = 0; i < x.size(); ++i) {
      hp[0][i] = 1 * std::exp(-0.5 * x[i] * x[i]); // Generate first two polynomilas manually
      hp[1][i] = 2 * x[i] * std::exp(-0.5 * x[i] * x[i]);
   }
   for (int n = 2; n < order; ++n) { // Then use recurrent chain
      for (size_t i = 0; i < x.size(); ++i) {
         hp[n][i] = (2 * x[i] * hp[n - 1][i] - 2 * (n - 1) * hp[n - 2][i]);
      }
   }
   return hp;
}

// Function to compute normalized physicists' Hermite polynomials for a given axis
std::vector<std::vector<float>> get_hermite(const HERMITE::OrderedVDF& data, int order, float vth,
                                            const std::vector<float>& u, int axis) {
   // Determine velocity limits for the given axis
   float v_min = data.v_limits[axis];
   float v_max = data.v_limits[axis + 3];
   // Generate velocity grid
   std::vector<float> v_axis = linspace(v_min, v_max, data.shape[axis]);
   for (auto& val : v_axis) {
      val = (val - u[axis]) / vth;
   }
   std::vector<std::vector<float>> hermite_vals = hermite(v_axis, order);
   for (int n = 0; n < order; ++n) {
      float norm_const = std::sqrt(std::pow(2, n) * factorial(n) * std::sqrt(M_PI) * vth);
      for (size_t i = 0; i < v_axis.size(); ++i) {
         hermite_vals[n][i] /= norm_const;
      }
   }
   return hermite_vals;
}

// get drift velocity from VDFdata structure: <-- Proper
std::vector<float> get_drift_velocity(const HERMITE::OrderedVDF& data) {
   std::vector<float> u(3);
   std::vector<float> vx = linspace(data.v_limits[0], data.v_limits[3], data.shape[0]);
   std::vector<float> vy = linspace(data.v_limits[1], data.v_limits[4], data.shape[1]);
   std::vector<float> vz = linspace(data.v_limits[2], data.v_limits[5], data.shape[2]);
   float n = 0;
   float dv = (data.v_limits[3] - data.v_limits[0]) / (data.shape[0]);
   for (size_t i = 0; i < data.shape[0]; ++i) {
      for (size_t j = 0; j < data.shape[1]; ++j) {
         for (size_t k = 0; k < data.shape[2]; ++k) {
            int index = i * data.shape[2] * data.shape[1] + j * data.shape[2] + k;
            u[0] += vx[i] * data.vdf_vals[index] * dv * dv * dv;
            u[1] += vy[j] * data.vdf_vals[index] * dv * dv * dv;
            u[2] += vz[k] * data.vdf_vals[index] * dv * dv * dv;
            n += data.vdf_vals[index] * dv * dv * dv;
         }
      }
   }
   for (float& val : u) { // Use reference to modify elements
      val /= n;
   }
   return u;
}

// get thermal velocity <-- Proper
float get_thermal_velocity(const HERMITE::OrderedVDF& data, std::vector<float> u) {
   float dv = (data.v_limits[3] - data.v_limits[0]) / (data.shape[0]);
   std::vector<float> vx = linspace(data.v_limits[0], data.v_limits[3], data.shape[0]);
   std::vector<float> vy = linspace(data.v_limits[1], data.v_limits[4], data.shape[1]);
   std::vector<float> vz = linspace(data.v_limits[2], data.v_limits[5], data.shape[2]);
   float n = 0;
   float Pxx = 0;
   float Pyy = 0;
   float Pzz = 0;
   for (size_t i = 0; i < data.shape[0]; ++i) {
      for (size_t j = 0; j < data.shape[1]; ++j) {
         for (size_t k = 0; k < data.shape[2]; ++k) {
            int index = i * data.shape[2] * data.shape[1] + j * data.shape[1] + k;
            Pxx += (vx[i] - u[0]) * ((vx[i] - u[0])) * data.vdf_vals[index] * dv * dv * dv;
            Pyy += (vy[j] - u[1]) * ((vy[j] - u[1])) * data.vdf_vals[index] * dv * dv * dv;
            Pzz += (vz[k] - u[2]) * ((vz[k] - u[2])) * data.vdf_vals[index] * dv * dv * dv;
            n += data.vdf_vals[index] * dv * dv * dv;
         }
      }
   }
   float vth = std::sqrt((Pxx + Pyy + Pzz) / (3 * n));
   return vth;
}

// calculate hermite spectra in 3D
std::vector<float> hermite_spectra_3d(const HERMITE::OrderedVDF& data, int order, float vth, const std::vector<float>& u) {
   std::vector<float> spectra(order * order * order);
   const float dv = (data.v_limits[3] - data.v_limits[0]) / (data.shape[0]);
   const std::vector<std::vector<float>> hermite_x = get_hermite(data, order, vth, u, 0);
   const std::vector<std::vector<float>> hermite_y = get_hermite(data, order, vth, u, 1);
   const std::vector<std::vector<float>> hermite_z = get_hermite(data, order, vth, u, 2);
   // loop over hermite
   for (int nx = 0; nx < order; ++nx) {
      for (int ny = 0; ny < order; ++ny) {
         for (int nz = 0; nz < order; ++nz) {
            const auto hermite_index = nx * (order) * (order) + ny * (order) + nz;
            float sum = 0.0f;
            for (size_t ix = 0; ix < data.shape[0]; ++ix) {
               for (size_t iy = 0; iy < data.shape[1]; ++iy) {
                  for (size_t iz = 0; iz < data.shape[0]; ++iz) {
                     const auto vspace_index=data.index(ix,iy,iz);
                     sum += data.vdf_vals[vspace_index] * hermite_x[nx][ix] * hermite_y[ny][iy] * hermite_z[nz][iz] *
                            dv * dv * dv;
                  }
               }
            }
            spectra[hermite_index] = sum;
         }
      }
   }
   return spectra;
}

// reconstruction
std::vector<float> reconstruct_vdf(const HERMITE::OrderedVDF& data, const std::vector<float>& spectra, int order, float vth,
                                   const std::vector<float>& u) {

   std::vector<float> f(data.shape[2] * data.shape[1] * data.shape[0]);
   const std::vector<std::vector<float>> hermite_x = get_hermite(data, order, vth, u, 0);
   const std::vector<std::vector<float>> hermite_y = get_hermite(data, order, vth, u, 1);
   const std::vector<std::vector<float>> hermite_z = get_hermite(data, order, vth, u, 2);
   const float dv = (data.v_limits[3] - data.v_limits[0]) / (data.shape[0]);
   float max_f = 0.0f;
   for (size_t vx = 0; vx < data.shape[0]; ++vx) {
      for (size_t vy = 0; vy < data.shape[1]; ++vy) {
         for (size_t vz = 0; vz < data.shape[2]; ++vz) {
            float sum = 0.0f;
            for (int nx = 0; nx < order; ++nx) {
               for (int ny = 0; ny < order; ++ny) {
                  for (int nz = 0; nz < order; ++nz) {
                     const std::size_t n = nx * (order) * (order) + ny * (order) + nz;
                     sum += spectra[n] * hermite_x[nx][vx] * hermite_y[ny][vy] * hermite_z[nz][vz];
                  }
               }
            }
            const std::size_t ind=data.index(vx,vy,vz);
            f[ind] = sum;
            if (sum > max_f)
               max_f = sum;
         }
      }
   }
   // // Apply thresholding
   // float threshold = 0.01f * max_f;
   // for (float& val : f) {
   //    if (val < threshold) val = 0.0f;
   // }

   return f;
}

HERMITE::HermSpectrum HERMITE::getHermiteSpectra(HERMITE::OrderedVDF vdfdata) {
   int order = 22; // define max order of the hermite decomposition
   std::vector<float> u = get_drift_velocity(vdfdata);
   float vth = get_thermal_velocity(vdfdata, u);
   std::vector<float> spectra = hermite_spectra_3d(vdfdata, order, vth, u);
   return HERMITE::HermSpectrum{.N_hermite_harmonic = order, .vth = vth, .u = u, .Spectrum = spectra};
}

HERMITE::OrderedVDF HERMITE::hermite_transform_back_and_forth(HERMITE::OrderedVDF vdfdata) {
   int order = 22; // define max order of the hermite decomposition
   // std::string fileName="vdf_41.bin"; // load vlsv vdf from binary
   // VDFdata data = read_vdf_bin(fileName);
   std::vector<float> u = get_drift_velocity(vdfdata);
   float vth = get_thermal_velocity(vdfdata, u);
   std::vector<float> spectra = hermite_spectra_3d(vdfdata, order, vth, u);
   std::vector<float> vdf_recon = reconstruct_vdf(vdfdata, spectra, order, vth, u);
   vdfdata.vdf_vals = vdf_recon;
   return vdfdata;
}

/// MAIN OVERWRITE
int HERMITE::overwrite_pop_spatial_cell_vdf(SpatialCell* sc, uint popID, const OrderedVDF& vdf) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   auto* blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer->size();
   const Real* blockParams = sc->get_block_parameters(popID);
   Realf* data = blockContainer->getData();
   Real sparse = getObjectWrapper().particleSpecies[popID].sparseMinValue;

   for (std::size_t n = 0; n < total_blocks; ++n) {
      auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      Realf* vdf_data = &data[n * WID3];
      for (uint k = 0; k < WID; ++k) {
         for (uint j = 0; j < WID; ++j) {
            for (uint i = 0; i < WID; ++i) {
               const Real dvx = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVX];
               const Real dvy = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVY];
               const Real dvz = (blockParams + BlockParams::N_VELOCITY_BLOCK_PARAMS)[BlockParams::DVZ];
               const std::size_t nx = std::ceil((vdf.v_limits[3] - vdf.v_limits[0]) / dvx);
               const std::size_t ny = std::ceil((vdf.v_limits[4] - vdf.v_limits[1]) / dvy);
               const std::size_t nz = std::ceil((vdf.v_limits[5] - vdf.v_limits[2]) / dvz);
               const Real vx = bp[BlockParams::VXCRD] + (i + 0.5) * bp[BlockParams::DVX];
               const Real vy = bp[BlockParams::VYCRD] + (j + 0.5) * bp[BlockParams::DVY];
               const Real vz = bp[BlockParams::VZCRD] + (k + 0.5) * bp[BlockParams::DVZ];
               const size_t bbox_i = std::min(static_cast<size_t>(std::floor((vx - vdf.v_limits[0]) / dvx)), nx - 1);
               const size_t bbox_j = std::min(static_cast<size_t>(std::floor((vy - vdf.v_limits[1]) / dvy)), ny - 1);
               const size_t bbox_k = std::min(static_cast<size_t>(std::floor((vz - vdf.v_limits[2]) / dvz)), nz - 1);
               const size_t index = bbox_i * (ny * nz) + bbox_j * nz + bbox_k;
               // vspace.at(index) += vdf_data[cellIndex(i, j, k)] / ratio;

               vdf_data[cellIndex(i, j, k)] = sparse * std::pow(10, vdf.vdf_vals.at(index));
               // vdf_data[cellIndex(i, j, k)] = vdf.vdf_vals.at(index);
            }
         }
      }
   } // over blocks
   return 0;
}


/// I need to extract and somehow write E and B
HERMITE::EBat HERMITE::dropB(spatial_cell::SpatialCell* sc){
	HERMITE::EBat eb;

      eb.B[0] = sc->parameters[CellParams::PERBXVOL] + sc->parameters[CellParams::BGBXVOL];
      eb.B[1] = sc->parameters[CellParams::PERBYVOL] + sc->parameters[CellParams::BGBYVOL];
      eb.B[2] = sc->parameters[CellParams::PERBZVOL] + sc->parameters[CellParams::BGBZVOL];

      eb.E[0] = sc->parameters[CellParams::EXVOL];
      eb.E[1] = sc->parameters[CellParams::EYVOL];
      eb.E[2] = sc->parameters[CellParams::EZVOL];

   return eb;
}




////// ROTATION of the VDF
// ROT LIMITS
std::array<double, 6> compute_rotated_limits(
    const std::vector<float>& vx,
    const std::vector<float>& vy,
    const std::vector<float>& vz,
    const Matrix3x3& R
) {
    float vmin[3] = {-3e6, -3e6, -3e6};
    float vmax[3] = {3e6, 3e6, 3e6};

    for (float x : vx) {
        for (float y : vy) {
            for (float z : vz) {
                Vector3 v = {x, y, z};
                Vector3 vr = {
                    R[0][0]*v[0] + R[0][1]*v[1] + R[0][2]*v[2],
                    R[1][0]*v[0] + R[1][1]*v[1] + R[1][2]*v[2],
                    R[2][0]*v[0] + R[2][1]*v[1] + R[2][2]*v[2]
                };
                for (int i = 0; i < 3; ++i) {
                    vmin[i] = std::min(vmin[i], vr[i]);
                    vmax[i] = std::max(vmax[i], vr[i]);
                }
            }
        }
    }

    return {vmin[0], vmin[1], vmin[2], vmax[0], vmax[1], vmax[2]};
}


// get ROTATED VDF structure for a given rotation matrix
HERMITE::OrderedVDF rotate_vdf(const HERMITE::OrderedVDF& input, const Matrix3x3& R) {
    HERMITE::OrderedVDF rotated = input;
    rotated.vdf_vals.assign(input.vdf_vals.size(), 0.0f);

    auto shape = input.shape;
    auto lin_x = linspace(input.v_limits[0], input.v_limits[3], shape[0]);
    auto lin_y = linspace(input.v_limits[1], input.v_limits[4], shape[1]);
    auto lin_z = linspace(input.v_limits[2], input.v_limits[5], shape[2]);

    for (size_t ix = 0; ix < shape[0]; ++ix) {
        for (size_t iy = 0; iy < shape[1]; ++iy) {
            for (size_t iz = 0; iz < shape[2]; ++iz) {
                Vector3 v_orig = {lin_x[ix], lin_y[iy], lin_z[iz]};
                Vector3 v_rot = {
                    R[0][0] * v_orig[0] + R[0][1] * v_orig[1] + R[0][2] * v_orig[2],
                    R[1][0] * v_orig[0] + R[1][1] * v_orig[1] + R[1][2] * v_orig[2],
                    R[2][0] * v_orig[0] + R[2][1] * v_orig[1] + R[2][2] * v_orig[2]
                };
                // Interpolate from input.vdf_vals at position v_rot
                float f_val = interpolate3D(input, v_rot); // You'll implement this
                size_t idx = rotated.index(ix, iy, iz);
                rotated.vdf_vals[idx] = f_val;
            }
        }
    }
    // Also rotate velocity limits!
    rotated.v_limits = compute_rotated_limits(lin_x, lin_y, lin_z, R);
    return rotated;
}


// get ROTATED VDF structure for a given spatial cell
//namespace HERMITE {
HERMITE::OrderedVDF HERMITE::rotate_vdf4cell(const HERMITE::OrderedVDF vdf_struct, spatial_cell::SpatialCell* sc) {
    HERMITE::OrderedVDF rotated_vdf = vdf_struct;
    rotated_vdf.vdf_vals.assign(vdf_struct.vdf_vals.size(), 0.0f);

        Vector3 E = {1.0f, 1.0f, 1.0f};
        Vector3 B = {1.0f, 2.0f, 1.0f};

        B[0] = sc->parameters[CellParams::PERBXVOL] + sc->parameters[CellParams::BGBXVOL];
        B[1] = sc->parameters[CellParams::PERBYVOL] + sc->parameters[CellParams::BGBYVOL];
        B[2] = sc->parameters[CellParams::PERBZVOL] + sc->parameters[CellParams::BGBZVOL];

        E[0] = sc->parameters[CellParams::EXVOL];
        E[1] = sc->parameters[CellParams::EYVOL];
        E[2] = sc->parameters[CellParams::EZVOL];

        Matrix3x3 R = build_rotation_matrix(E, B);
//	std::cout << "Rotation matrix R:\n";
//	for (int i=0; i<3; i++) {
//	    std::cout << R[i][0] << " " << R[i][1] << " " << R[i][2] << "\n";
//	}

//    LIMITS for the old implementation of rotation	
//    std::vector<float>  lin_x = linspace(vdf_struct.v_limits[0], vdf_struct.v_limits[3], vdf_struct.shape[0]);
//    std::vector<float>  lin_y = linspace(vdf_struct.v_limits[1], vdf_struct.v_limits[4], vdf_struct.shape[1]);
//    std::vector<float>  lin_z = linspace(vdf_struct.v_limits[2], vdf_struct.v_limits[5], vdf_struct.shape[2]);
    // 1. Calculate new limits and axes
    float max_abs_v = std::max({
        std::abs(vdf_struct.v_limits[0]), std::abs(vdf_struct.v_limits[3]),
        std::abs(vdf_struct.v_limits[1]), std::abs(vdf_struct.v_limits[4]),
        std::abs(vdf_struct.v_limits[2]), std::abs(vdf_struct.v_limits[5])
    });
    float scale_factor = 1.3f;
    float new_limit = max_abs_v * scale_factor;

    std::array<float, 6> new_v_limits = {
        -new_limit, -new_limit, -new_limit,
         new_limit,  new_limit,  new_limit
    };

    auto nx = vdf_struct.shape[0];
    auto ny = vdf_struct.shape[1];
    auto nz = vdf_struct.shape[2];

    auto lin_x = linspace(new_v_limits[0], new_v_limits[3], nx);
    auto lin_y = linspace(new_v_limits[1], new_v_limits[4], ny);
    auto lin_z = linspace(new_v_limits[2], new_v_limits[5], nz);

//	std::cout << "v_limits: ";
//	for (auto v : vdf_struct.v_limits) std::cout << v << " ";
//	std::cout << "\n";

//	std::cout << "lin_x size: " << lin_x.size() << " values: ";
//	for (auto x : lin_x) std::cout << x << " ";
//	std::cout << "\n";

    // 2. Compute inverse rotation matrix (transpose of R for rotation matrices)
    std::array<std::array<float, 3>, 3> R_inv;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        R_inv[i][j] = R[j][i];

    // 3. For each point in new grid:
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                // Coordinate in rotated frame
                float vx_rot = lin_x[i];
                float vy_rot = lin_y[j];
                float vz_rot = lin_z[k];

                // Apply inverse rotation to get original frame velocity
                float vx_orig = R_inv[0][0] * vx_rot + R_inv[0][1] * vy_rot + R_inv[0][2] * vz_rot;
                float vy_orig = R_inv[1][0] * vx_rot + R_inv[1][1] * vy_rot + R_inv[1][2] * vz_rot;
                float vz_orig = R_inv[2][0] * vx_rot + R_inv[2][1] * vy_rot + R_inv[2][2] * vz_rot;

                // Interpolate in original VDF
                rotated_vdf.vdf_vals[idx(rotated_vdf, i, j, k)] = trilinear_interpolate(vdf_struct, vx_orig, vy_orig, vz_orig);
            }
        }
    }    

    return rotated_vdf;
}

//} // end Hermite namespace

HERMITE::HermSpectrum HERMITE::getHERMITE_VDFRot(spatial_cell::SpatialCell* sc, HERMITE::OrderedVDF original_vdf ){

  	HERMITE::HermSpectrum SpectrStruct;
        Vector3 E = {1.0f, 1.0f, 1.0f};
        Vector3 B = {1.0f, 2.0f, 1.0f};

        B[0] = sc->parameters[CellParams::PERBXVOL] + sc->parameters[CellParams::BGBXVOL];
        B[1] = sc->parameters[CellParams::PERBYVOL] + sc->parameters[CellParams::BGBYVOL];
        B[2] = sc->parameters[CellParams::PERBZVOL] + sc->parameters[CellParams::BGBZVOL];

        E[0] = sc->parameters[CellParams::EXVOL];
        E[1] = sc->parameters[CellParams::EYVOL];
        E[2] = sc->parameters[CellParams::EZVOL];

        Matrix3x3 R = build_rotation_matrix(E, B);        
	auto rotated_vdf = rotate_vdf(original_vdf, R);

	int order = 22;
	auto u_rot = get_drift_velocity(rotated_vdf);
	auto vth_rot = get_thermal_velocity(rotated_vdf, u_rot);
	std::vector<float> spectrum_rot = hermite_spectra_3d(rotated_vdf, order, vth_rot, u_rot);	

	return HERMITE::HermSpectrum{.N_hermite_harmonic = order, .vth = vth_rot, .u = u_rot, .Spectrum = spectrum_rot};	
}





