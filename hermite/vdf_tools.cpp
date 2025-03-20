#include "../object_wrapper.h"
#include "../spatial_cell_wrapper.hpp"
#include "../velocity_blocks.h"
#include <array>
#include "vdf_tools.h"




HERMITE::OrderedVDF  HERMITE::extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed( SpatialCell* sc, uint popID,
                                                                                       int zoom) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   if (zoom != 1) {
      throw std::runtime_error("Zoom is not supported yet!");
   }
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer.size();
   const Real* blockParams = sc->get_block_parameters(popID);

   Real sparse = getObjectWrapper().particleSpecies[popID].sparseMinValue;

   // xmin,ymin,zmin,xmax,ymax,zmax;
   std::array<Real, 6> vlims{std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::max(),
                             std::numeric_limits<Real>::max(),    std::numeric_limits<Real>::lowest(),
                             std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest()};

   // This pass is computing the active vmesh limits
   // Store dvx,dvy,dvz here
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

   assert(isPow2(static_cast<size_t>(std::abs(zoom))));
   float ratio = (zoom > 0) ? static_cast<float>(std::abs(zoom)) : 1.0 / static_cast<float>(std::abs(zoom));
   assert(ratio > 0);

   const Real target_dvx = dvx * ratio;
   const Real target_dvy = dvy * ratio;
   const Real target_dvz = dvz * ratio;
   std::size_t nx = std::ceil((vlims[3] - vlims[0]) / target_dvx);
   std::size_t ny = std::ceil((vlims[4] - vlims[1]) / target_dvy);
   std::size_t nz = std::ceil((vlims[5] - vlims[2]) / target_dvz);
   // printf("VDF min box is %zu , %zu %zu \n ", nx, ny, nz);
   
   std::unordered_set<vmesh::GlobalID> ignore_list;
   for (std::size_t k = 0; k < nz; ++k) {
      for (std::size_t j = 0; j < ny; ++j) {
         for (std::size_t i = 0; i < nx; ++i) {
               const Real vx = vlims[0] + (i + 0.5) * dvx;
               const Real vy = vlims[1] + (j + 0.5) * dvy;
               const Real vz = vlims[2] + (k + 0.5) * dvz;
               const std::array<Real,3>coords={vx,vy,vz};
               const auto gid=sc->get_velocity_block(popID, &coords[0]);
               ignore_list.insert(gid);
         }
      }
   }
   
   Realf* data = blockContainer.getData();
   std::vector<Realf> vspace(nx * ny * nz, Realf(0));
   for (std::size_t n = 0; n < total_blocks; ++n) {
      const auto bp = blockParams + n * BlockParams::N_VELOCITY_BLOCK_PARAMS;
      const Realf* vdf_data = &data[n * WID3];
      const vmesh::GlobalID gid = sc->get_velocity_block_global_id(n, popID);
      (void)ignore_list.erase(gid);
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
                  if (vdf_data[cellIndex(i, j, k)]<sparse){
                        vspace.at(index)=0.0;
                     } else {                        
                        vspace.at(index) += vdf_data[cellIndex(i, j, k)] / ratio;                         
                     }                     
                  //vspace.at(index) += vdf_data[cellIndex(i, j, k)] / ratio; // that my sh**
               } else {
                  // Same value in all bins
                  int max_off = 1 / ratio;
                  for (int off_z = 0; off_z <= max_off; off_z++) {
                     for (int off_y = 0; off_y <= max_off; off_y++) {
                        for (int off_x = 0; off_x <= max_off; off_x++) {
                           const size_t index = (bbox_i + off_x) * (ny * nz) + (bbox_j + off_y) * nz + (bbox_k + off_z);
                           if (index < vspace.size()) {
                              if (vdf_data[cellIndex(i, j, k)]<sparse){
                                 vspace.at(index)=0.0;
                                } else {
                                 vspace.at(index) = vdf_data[cellIndex(i, j, k)]; //it was only this thing here without sparcity
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

   std::vector<vmesh::GlobalID >ignored;
   ignored.reserve(ignore_list.size());
   for (auto it = ignore_list.begin(); it != ignore_list.end(); ) {
      ignored.push_back(std::move(ignore_list.extract(it++).value()));
   }
 
   for (uint i=0; i<vspace.size(); ++i){
       vspace[i] = std::log10(std::max(vspace[i], static_cast<float>(sparse))) - std::log10(sparse);
    }

   return HERMITE::OrderedVDF{.blocks_to_ignore=ignored,.sparse_vdf_bytes=total_blocks*WID*WID*WID*sizeof(Realf),.vdf_vals = vspace, .v_limits = vlims, .shape = {nx, ny, nz}};
}


void dump_vdf_to_binary_file(const char* filename,uint popID, CellID cid,
   dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) {
   SpatialCell* sc = mpiGrid[cid];
   assert(sc && "Invalid Pointer to Spatial Cell !");
   HERMITE::OrderedVDF vdf = HERMITE::extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed(sc, popID, 1);
   // vdf.save_to_file(filename);
}


//  Hermite decomposition functions are below 
// factorial
unsigned long long factorial(int n) {
   if (n < 0) {
       throw std::invalid_argument("Factorial is not defined for negative numbers.");
   }
   unsigned long long result = 1;
   for (int i = 2; i <= n; ++i) {
       result *= i;       
       if (result < i) {
           throw std::overflow_error("Factorial result is too large!");
       }
   }
   return result;
}

// linspace
std::vector<float> linspace(double start, double end, int len){
 float step = (end-start)/(len - 1);
 std::vector<float> x(len);
 for(int i=0; i<len; ++i){
   x[i] = start + i*step;
   }
 return x;
}

// generate hermite polinomials up to the given order
std::vector<std::vector<float>> hermite(std::vector<float>& x, int order){
   // Recurrence relation: H_n+1 = 2*x*H_n -2*n*H_n-1
   // H0 = 1, H1 = 2x
   // base function with Gauss weights: H_n * exp(-0.5*v^2)
     std::vector<std::vector<float>> hp(order, std::vector<float>(x.size()) ) ;
     for(size_t i=0; i<x.size(); ++i){
       hp[0][i] = 1 * std::exp(-0.5*x[i]*x[i]) ;       //Generate first two polynomilas manually
       hp[1][i] = 2*x[i] * std::exp(-0.5*x[i]*x[i]);
     }
     for(int n=2; n<order; ++n){  // Then use recurrent chain
       for(size_t i=0; i<x.size(); ++i){     
         hp[n][i] = (2*x[i]*hp[n-1][i] - 2*(n-1)*hp[n-2][i] ) ;
       }
     }
     return hp;
   }


   // Function to compute normalized physicists' Hermite polynomials for a given axis
   std::vector<std::vector<float>> get_hermite(const HERMITE::OrderedVDF& data, int order, float vth, const std::vector<float>& u, int axis) {
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
   
 


//get drift velocity from VDFdata structure: <-- Proper
std::vector<float> get_drift_velocity(HERMITE::OrderedVDF data){   
   std::vector<float> u(3);
   std::vector<float> vx = linspace(data.v_limits[0],data.v_limits[3],data.shape[0]);
   std::vector<float> vy = linspace(data.v_limits[1],data.v_limits[4],data.shape[1]);
   std::vector<float> vz = linspace(data.v_limits[2],data.v_limits[5],data.shape[2]);
   float n = 0;
   float dv = (data.v_limits[3] - data.v_limits[0]) / (data.shape[0] );
   for(size_t i=0; i< data.shape[0]; ++i){
      for(size_t j=0; j< data.shape[1]; ++j){
        for(size_t k=0; k< data.shape[2]; ++k){          
          int index = i*data.shape[2]*data.shape[1] + j*data.shape[2] + k;
          u[0] +=  vx[i] * data.vdf_vals[index] * dv * dv * dv;
          u[1] +=  vy[j] * data.vdf_vals[index] * dv * dv * dv;
          u[2] +=  vz[k] * data.vdf_vals[index] * dv * dv * dv;
          n += data.vdf_vals[index] * dv * dv * dv;
         }
      }
   }
   for (float& val : u) {  // Use reference to modify elements
      val /= n;
  }
   return u;
 }
 
 //get thermal velocity <-- Proper 
 float get_thermal_velocity(HERMITE::OrderedVDF data, std::vector<float> u){
   float dv = (data.v_limits[3] - data.v_limits[0]) / (data.shape[0] );
   std::vector<float> vx = linspace(data.v_limits[0],data.v_limits[3],data.shape[0]);
   std::vector<float> vy = linspace(data.v_limits[1],data.v_limits[4],data.shape[1]);
   std::vector<float> vz = linspace(data.v_limits[2],data.v_limits[5],data.shape[2]);
   float n = 0;
   float Pxx = 0;
   float Pyy = 0;
   float Pzz = 0;
   for(size_t i=0; i<data.shape[0]; ++i){
      for(size_t j=0; j<data.shape[1]; ++j){
        for(size_t k=0; k<data.shape[2]; ++k){
         int index = i*data.shape[2]*data.shape[1] + j*data.shape[1] + k;
         Pxx += (vx[i] - u[0]) * ((vx[i] - u[0])) * data.vdf_vals[index] * dv * dv * dv ;
         Pyy += (vy[j] - u[1]) * ((vy[j] - u[1])) * data.vdf_vals[index] * dv * dv * dv ;
         Pzz += (vz[k] - u[2]) * ((vz[k] - u[2])) * data.vdf_vals[index] * dv * dv * dv ;
         n += data.vdf_vals[index] * dv * dv * dv;
        }
      }
   }
  float vth = std::sqrt(  (Pxx+Pyy+Pzz) / (3 * n) );
return vth;
}

 
// calculate hermite spectra in 3D
std::vector<float> hermite_spectra_3d(HERMITE::OrderedVDF data, int order, float vth, std::vector<float> u){
   std::vector<float> spectra(order*order*order);
   int hermite_index;
   int vspace_index;
   float sum;
   float dv = (data.v_limits[3]-data.v_limits[0]) / (data.shape[0])  ;
   std::vector<std::vector<float>> hermite_x = get_hermite(data, order, vth, u, 0);
   std::vector<std::vector<float>> hermite_y = get_hermite(data, order, vth, u, 1);
   std::vector<std::vector<float>> hermite_z = get_hermite(data, order, vth, u, 2);
   // loop over hermite
   for(int nx=0; nx<order; ++nx){
     for(int ny=0; ny<order; ++ny){
       for(int nz=0; nz<order; ++nz){
         hermite_index = nx*(order)*(order)+ny*(order)+nz;
         //loop over vspace   
         sum=0;
         for(size_t ix=0; ix<data.shape[0]; ++ix){
           for(size_t iy=0; iy<data.shape[1]; ++iy){
             for(size_t iz=0; iz<data.shape[0]; ++iz){
               vspace_index=ix*(data.shape[2])*(data.shape[1])+iy*(data.shape[2])+iz;
               sum+=data.vdf_vals[vspace_index]*hermite_x[nx][ix]*hermite_y[ny][iy]*hermite_z[nz][iz]*dv*dv*dv;
             }
           }
         }
         spectra[hermite_index]=sum;
       }
     }
   }
   return spectra;
 }

// reconstruction
std::vector<float> reconstruct_vdf(HERMITE::OrderedVDF data, std::vector<float> spectra, int order, float vth, std::vector<float> u){

   std::vector<float> f(data.shape[2]*data.shape[1]*data.shape[0]);
   
   std::vector<std::vector<float>> hermite_x = get_hermite(data, order, vth, u, 0) ; 
   std::vector<std::vector<float>> hermite_y = get_hermite(data, order, vth, u, 1) ;
   std::vector<std::vector<float>> hermite_z = get_hermite(data, order, vth, u, 2) ;
   
   float dv = (data.v_limits[3]-data.v_limits[0]) / (data.shape[0]);
   float max_f = 0.0f; 

   for (size_t vx=0; vx<data.shape[0]; ++vx){    
     for (size_t vy=0; vy<data.shape[1]; ++vy){
       for (size_t vz=0; vz<data.shape[2]; ++vz){
         int ind = vx*data.shape[1]*data.shape[2] + vy*data.shape[2] + vz;
         float sum = 0.0f;
         for(int nx=0; nx<order; ++nx){
           for(int ny=0; ny<order; ++ny){
             for(int nz=0; nz<order; ++nz){
               int n = nx*(order)*(order) + ny*(order) + nz;
               sum += spectra[n]*hermite_x[nx][vx]*hermite_y[ny][vy]*hermite_z[nz][vz];
             }
           }
         }
         f[ind]=sum;
         if (sum > max_f) max_f = sum;
       }
     }
   }
   // Apply thresholding 
   float threshold = 0.05f * max_f;
   for (float& val : f) {
      if (val < threshold) val = 0.0f;
   }

 return f;
 }

 
 HERMITE::HermSpectrum HERMITE::getHermiteSpectra(HERMITE::OrderedVDF vdfdata){
   int order=18; // define max order of the hermite decomposition
   std::vector<float> u = get_drift_velocity(vdfdata); 
   float vth = get_thermal_velocity(vdfdata, u);
   std::vector<float> spectra=hermite_spectra_3d(vdfdata, order, vth, u); 
   return HERMITE::HermSpectrum{.N_hermite_harmonic = order, .vth = vth, .u = u, .HermSpectrum = spectra };
 }


 HERMITE::OrderedVDF HERMITE::hermite_transform_back_and_forth(HERMITE::OrderedVDF vdfdata){
   int order=18; // define max order of the hermite decomposition
   // std::string fileName="vdf_41.bin"; // load vlsv vdf from binary
   // VDFdata data = read_vdf_bin(fileName);
   std::vector<float> u = get_drift_velocity(vdfdata); 
   float vth = get_thermal_velocity(vdfdata, u);
   std::vector<float> spectra=hermite_spectra_3d(vdfdata, order, vth, u);  
   std::vector<float> vdf_recon = reconstruct_vdf( vdfdata, spectra, order, vth, u);
   vdfdata.vdf_vals = vdf_recon;
   return vdfdata;
}

/// MAIN OVERWRITE
int HERMITE::overwrite_pop_spatial_cell_vdf(SpatialCell* sc, uint popID, const OrderedVDF& vdf) {
   assert(sc && "Invalid Pointer to Spatial Cell !");
   vmesh::VelocityBlockContainer<vmesh::LocalID>& blockContainer = sc->get_velocity_blocks(popID);
   const size_t total_blocks = blockContainer.size();
   const Real* blockParams = sc->get_block_parameters(popID);
   Realf* data = blockContainer.getData();
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

               vdf_data[cellIndex(i, j, k)] = sparse * std::pow(10,vdf.vdf_vals.at(index));
               
            }
         }
      }
   } // over blocks
   return 0;
}

