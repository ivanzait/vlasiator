/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "Harris.h"

using namespace spatial_cell;

namespace projects {
Harris::Harris() : TriAxisSearch() {}
Harris::~Harris() {}

bool Harris::initialize(void) { return Project::initialize(); }

void Harris::addParameters() {
   typedef Readparameters RP;
   RP::add("Harris.currentSheetType", "Type of the current sheet initialization: 1. single layer, 2. double layer", 1);
   RP::add("Harris.scale", "Harris sheet scale size (m)", 150000.0);
   RP::add("Harris.Bx0", "Reference Magnetic field (T)", 8.33061003094e-8);
   RP::add("Harris.By0", "Reference Magnetic field (T)", 8.33061003094e-8);
   RP::add("Harris.Bz0", "Reference Magnetic field (T)", 8.33061003094e-8);
   RP::add("Harris.Psi0", "Perturbation of the magnetic flux (T*m)", 0.0);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;

      RP::add(pop + "_Harris.n", "Reference number density (m^-3)", 1.0e6);
      RP::add(pop + "_Harris.T", "Species temperature", 1.0);
      RP::add(pop + "_Harris.nSpaceSamples", "Number of sampling points per spatial dimension.", 2);
      RP::add(pop + "_Harris.nVelocitySamples", "Number of sampling points per velocity dimension.", 2);
   }
}

void Harris::getParameters() {
   Project::getParameters();
   typedef Readparameters RP;
   RP::get("Harris.currentSheetType", currentSheetType);
   RP::get("Harris.scale", lambda);
   RP::get("Harris.Bx0", Bx0);
   RP::get("Harris.By0", By0);
   RP::get("Harris.Bz0", Bz0);
   RP::get("Harris.Psi0", Psi0);

   if (currentSheetType != 1 && currentSheetType != 2) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         std::cerr << "unknown initialization type: " << currentSheetType << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      HarrisSpeciesParameters sP;

      RP::get(pop + "_Harris.n", sP.n);
      RP::get(pop + "_Harris.T", sP.T);
      RP::get(pop + "_Harris.nSpaceSamples", sP.nSpaceSamples);
      RP::get(pop + "_Harris.nVelocitySamples", sP.nVelocitySamples);

      speciesParams.push_back(sP);
   }
}

Real Harris::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                           creal& dvz, const uint popID) const {
   const HarrisSpeciesParameters& sP = speciesParams[popID];
   Real mass = getObjectWrapper().particleSpecies[popID].mass;
   auto KB = physicalconstants::K_B;

   Real f, n;

   if (currentSheetType == 1) {
      Real sech2 = sqr(1.0 / std::cosh(z / lambda));
      n = sP.n * sech2 + 0.2 * sP.n;
   } else if (currentSheetType == 2) {
      const Real Lx = P::xmax - P::xmin;
      const Real Lz = P::zmax - P::zmin;

      n = sP.n * (sqr(1.0 / std::cosh((z - 0.25 * Lz) / lambda)) + sqr(1.0 / std::cosh((z + 0.25 * Lz) / lambda))) +
          0.2 * sP.n;
   }

   f = n * std::pow(mass / (2.0 * M_PI * KB * sP.T), 1.5) *
       exp(-mass * (sqr(vx) + sqr(vy) + sqr(vz)) / (2.0 * KB * sP.T));

   return f;
}


Real Harris::getBiMaxwellian(const uint popID, creal rho, creal Tpar, creal Tperp, creal vpar, creal vperp) const {
   const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
   Real c1 = MASS / (2.0 * physicalconstants::K_B);
   Real c2 = c1 / M_PI;
   Real f = rho / (sqrt(Tpar) * Tperp) * sqrt(c2) * c2 * exp(-c1 * (vperp * vperp) / Tperp - c1 * vpar * vpar / Tpar);
   return f;
}




Realf Harris::fillPhaseSpace(spatial_cell::SpatialCell *cell,
                                       const uint popID,
                                       const uint nRequested
      ) const {
      const HarrisSpeciesParameters& sP = speciesParams[popID];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
//      Real initRho = this->getCorrectNumberDensity(cell, popID);
	       	 
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];
      Real initRho = sP.n * (0.2 + 1.0 / pow(cosh(z / (lambda)), 2.0) )  ;
      Real initT = sP.T;
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

      #ifdef USE_GPU
      vmesh::VelocityMesh *vmesh = cell->dev_get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
      #else
      vmesh::VelocityMesh *vmesh = cell->get_velocity_mesh(popID);
      vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
      #endif

//      if (emptyBox == true) {
//         Realf* bufferData = cell->get_velocity_blocks(popID)->getData();
  //       std::memset(bufferData, 0, nRequested*WID3*sizeof(Realf));
    //     return 0;
    //  }   

      // Loop over blocks
      Realf rhosum = 0;
      arch::parallel_reduce<arch::null>(
         {WID, WID, WID, nRequested},
         ARCH_LOOP_LAMBDA (const uint i, const uint j, const uint k, const uint initIndex, Realf *lsum ) { 
            vmesh::GlobalID *GIDlist = vmesh->getGrid()->data();
            Realf* bufferData = VBC->getData();
            const vmesh::GlobalID blockGID = GIDlist[initIndex];
            // Calculate parameters for new block
            Real blockCoords[6];
            vmesh->getBlockInfo(blockGID,&blockCoords[0]);
            creal vxBlock = blockCoords[0];
            creal vyBlock = blockCoords[1];
            creal vzBlock = blockCoords[2];
            creal dvxCell = blockCoords[3];
            creal dvyCell = blockCoords[4];
            creal dvzCell = blockCoords[5];
            ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
               creal vx = vxBlock + (i+0.5)*dvxCell - initV0X;
               creal vy = vyBlock + (j+0.5)*dvyCell - initV0Y;
               creal vz = vzBlock + (k+0.5)*dvzCell - initV0Z;
               const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
               bufferData[initIndex*WID3 + k*WID2 + j*WID + i] = value;
               //lsum[0] += value;
            };  
         }, rhosum);
      return rhosum;
   } 

  /* Evaluates local SpatialCell properties for the project and population,
      then evaluates the phase-space density at the given coordinates.
      Used as a probe for projectTriAxisSearch.
   */
   Realf Harris::probePhaseSpace(spatial_cell::SpatialCell *cell,
                                        const uint popID,
                                        Real vx_in, Real vy_in, Real vz_in
      ) const {
      const HarrisSpeciesParameters& sP = speciesParams[popID];
      // Fetch spatial cell center coordinates
      const Real x  = cell->parameters[CellParams::XCRD] + 0.5*cell->parameters[CellParams::DX];
      const Real y  = cell->parameters[CellParams::YCRD] + 0.5*cell->parameters[CellParams::DY];
      const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];

      const Real mass = getObjectWrapper().particleSpecies[popID].mass;
//      Real initRho = sP.n;

//       const Real z  = cell->parameters[CellParams::ZCRD] + 0.5*cell->parameters[CellParams::DZ];
       Real initRho = sP.n * (0.2 + 1.0 / pow(cosh(z / (lambda)), 2.0) )  ;

      Real initT = sP.T;
      // Note: bulk V is zero, according to this and getV0().
      const Real initV0X = 0;
      const Real initV0Y = 0;
      const Real initV0Z = 0;

   //   initRho *= mass * (1.0 + 5.0 / pow(cosh(z / (lambda)), 2.0));
      creal vx = vx_in - initV0X;
      creal vy = vy_in - initV0Y;
      creal vz = vz_in - initV0Z;
      const Realf value = MaxwellianPhaseSpaceDensity(vx,vy,vz,initT,initRho,mass);
      return value;
   }



void Harris::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

std::vector<std::array<Real, 3>> Harris::getV0(creal x, creal y, creal z, const uint popID) const {
   std::vector<std::array<Real, 3>> V0;
   std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
   V0.push_back(v);
   return V0;
}


void Harris::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                              FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                              FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   setBackgroundFieldToZero(BgBGrid);

   const Real Lx = P::xmax - P::xmin;
   const Real Lz = P::zmax - P::zmin;

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

      if (currentSheetType == 1) {
#pragma omp parallel for collapse(3)
         for (int i = 0; i < localSize[0]; ++i) {
            for (int j = 0; j < localSize[1]; ++j) {
               for (int k = 0; k < localSize[2]; ++k) {
                  const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);

                  cell->at(fsgrids::bfield::PERBX) = Bx0 * tanh((x[2] + 0.5 * perBGrid.DZ) / lambda) +
                                                     Psi0 * (-M_PI / Lz) *
                                                         cos(2 * M_PI * (x[0] + 0.5 * perBGrid.DX) / Lx) *
                                                         sin(M_PI * (x[2] + 0.5 * perBGrid.DZ) / Lz);
                  cell->at(fsgrids::bfield::PERBY) = 0.0;
                  cell->at(fsgrids::bfield::PERBZ) = Psi0 * (-2 * M_PI / Lx) *
                                                     sin(2 * M_PI * (x[0] + 0.5 * perBGrid.DX) / Lx) *
                                                     cos(M_PI * (x[2] + 0.5 * perBGrid.DZ) / Lz);
               }
            }
         }
      } else if (currentSheetType == 2) {
#pragma omp parallel for collapse(3)
         for (int i = 0; i < localSize[0]; ++i) {
            for (int j = 0; j < localSize[1]; ++j) {
               for (int k = 0; k < localSize[2]; ++k) {
                  const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);

                  cell->at(fsgrids::bfield::PERBX) =
                      Bx0 * (tanh((x[2] + 0.5 * perBGrid.DZ - 0.25 * Lz) / lambda) -
                             tanh((x[2] + 0.5 * perBGrid.DZ + 0.25 * Lz) / lambda) + 1.0) +
                      Psi0 * (-M_PI / Lz) * cos(2 * M_PI * (x[0] + 0.5 * perBGrid.DX) / Lx) *
                          sin(2 * M_PI * (x[2] + 0.5 * perBGrid.DZ) / Lz - 0.25);
                  cell->at(fsgrids::bfield::PERBY) = 0.0;
                  cell->at(fsgrids::bfield::PERBZ) = Psi0 * (-2 * M_PI / Lx) *
                                                     sin(2 * M_PI * (x[0] + 0.5 * perBGrid.DX) / Lx) *
                                                     cos(2 * M_PI * (x[2] + 0.5 * perBGrid.DZ) / Lz);
               }
            }
         }
      }
   }
}

} // namespace projects
