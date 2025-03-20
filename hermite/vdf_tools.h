#pragma once
#include "stdlib.h"
#include <algorithm>
#include <array>
#include <vector>
#include <cstdint>
#include <fstream>
#include "../definitions.h"
#include <iostream>
#include "../spatial_cell_cpu.hpp"

namespace HERMITE {

   struct HermSpectrum {
      int N_hermite_harmonic;
      float vth;
      std::vector<float> u;
      std::vector<float> HermSpectrum;

      // std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
      //    return i * (N_hermite_harmonic * N_hermite_harmonic) + j * N_hermite_harmonic + k;
      // }
      // Realf& at(std::size_t i, std::size_t j, std::size_t k) noexcept { return HermSpectrum.at(index(i, j, k)); }   
      // const Realf& at(std::size_t i, std::size_t j, std::size_t k) const noexcept { return HermSpectrum.at(index(i, j, k)); }

      bool save_to_file(const char* filename) const noexcept {            
         std::ofstream file(filename, std::ios::binary);
         if (!file) {
            std::cerr << "Could not open file for writting! Exiting!" << std::endl;
            return false;
         }
         file.write( reinterpret_cast<const char*>(&N_hermite_harmonic), sizeof(int));         
         if (!file) {
            std::cerr << "Error writing Number of harmonics!" << std::endl;
            return false;
         }
         file.write( reinterpret_cast<const char*>(&vth), sizeof(float) );         
         if (!file) {
            std::cerr << "Error writing thermal speed !" << std::endl;
            return false;
         }
         file.write( reinterpret_cast<const char*>(u.data()), 3*sizeof(float));         
         if (!file) {
            std::cerr << "Error writing bulk velocity !" << std::endl;
            return false;
         }
         file.write( reinterpret_cast<const char*>(HermSpectrum.data()), pow(N_hermite_harmonic,3)*sizeof(float));         
         if (!file) {
            std::cerr << "Error writing vector of Hermite harmonics !" << std::endl;
            return false;
         }
      return true;
      } // end saving spectrum data

   }; // end HermSpectrum structure
   

   struct OrderedVDF {
      
      std::vector<vmesh::GlobalID> blocks_to_ignore;
      std::size_t sparse_vdf_bytes = {0};
      std::vector<Realf> vdf_vals;
      std::array<Real, 6> v_limits;     // vx_min,vy_min,vz_min,vx_max,vy_max,vz_max
      std::array<std::size_t, 3> shape; // x,y,z
      std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
         return i * (shape[1] * shape[2]) + j * shape[2] + k;
      }
   
      Realf& at(std::size_t i, std::size_t j, std::size_t k) noexcept { return vdf_vals.at(index(i, j, k)); }
   
      const Realf& at(std::size_t i, std::size_t j, std::size_t k) const noexcept { return vdf_vals.at(index(i, j, k)); }
   
      bool save_to_file(const char* filename) const noexcept {
         std::ofstream file(filename, std::ios::out | std::ios::binary);
         if (!file) {
            std::cerr << "Could not open file for writting! Exiting!" << std::endl;
            return false;
         }
         file.write((char*)shape.data(), 3 * sizeof(size_t));
         if (!file) {
            std::cerr << "Error writing shape data to file!" << std::endl;
            return false;
         }

         file.write((char*)v_limits.data(), 6 * sizeof(Real));
         if (!file) {
            std::cerr << "Error writing shape data to file!" << std::endl;
            return false;
         }         

         file.write((char*)vdf_vals.data(), vdf_vals.size() * sizeof(Realf));
         if (!file) {
            std::cerr << "Error writing vdf_vals data to file!" << std::endl;
            return false;
         }
         return true;
      }
   };


   

////// V COOORDINATES - DO I NEED ???

   struct VCoords {
      Real vx, vy, vz;
      VCoords operator+(const VCoords& other) { return {vx + other.vx, vy + other.vy, vz + other.vz}; }
      VCoords operator-(const VCoords& other) { return {vx - other.vx, vy - other.vy, vz - other.vz}; }
   };




OrderedVDF extract_pop_vdf_from_spatial_cell_ordered_min_bbox_zoomed( SpatialCell* sc, uint popID,int zoom);

HermSpectrum getHermiteSpectra(OrderedVDF vdfdata);

OrderedVDF hermite_transform_back_and_forth(OrderedVDF vdfdata);

 int overwrite_pop_spatial_cell_vdf(SpatialCell* sc, uint popID, const OrderedVDF& vdf);

// int overwrite_pop_spatial_cell_vdf(SpatialCell* sc, uint popID, const std::vector<Realf>& new_vspace);


} // End Hermite namespace 



