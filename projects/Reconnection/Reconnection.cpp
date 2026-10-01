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

/*
 * 2D GEM-like magnetic reconnection project (Birn et al. 2001).
 * See Reconnection.h for geometry / physics description.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "Reconnection.h"

using namespace spatial_cell;

namespace projects {

Reconnection::Reconnection() : TriAxisSearch() {}
Reconnection::~Reconnection() {}

bool Reconnection::initialize(void) { return Project::initialize(); }

// ---------------------------------------------------------------------------
void Reconnection::addParameters() {
   typedef Readparameters RP;
   RP::add<Real>("Reconnection.lambda",
           "Current-sheet half-thickness (m)",
           this->lambda,
           226284.5); // ~ 1 d_i for n=1e7, B=10 nT
   RP::add<Real>("Reconnection.Bx0",
           "Lobe (asymptotic) magnetic field magnitude Bx (T)",
           this->Bx0,
           1.0e-8);
   RP::add<Real>("Reconnection.By0",
           "Guide field By (T); 0 = anti-parallel reconnection",
           this->By0,
           0.0);
   RP::add<Real>("Reconnection.Bz0",
           "Initial normal field Bz (T); usually 0",
           this->Bz0,
           0.0);
   RP::add<Real>("Reconnection.Psi0",
           "Perturbation amplitude of the magnetic flux function (T*m); "
           "psi = Psi0*cos(2pi*x/Lx)*cos(pi*z/Lz)",
           this->Psi0,
           0.0);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      ReconnectionSpeciesParameters* sP = new ReconnectionSpeciesParameters();
      speciesParamsRead.push_back(sP);

      RP::add<Real>(pop + "_Reconnection.n",
              "Peak number density (m^-3)", sP->n, 1.0e7);
      RP::add<Real>(pop + "_Reconnection.T",
              "Temperature (K)", sP->T, 2.0e6);
      RP::add<Real>(pop + "_Reconnection.n_bg_frac",
              "Background density fraction n_bg/n0 (dimensionless)", sP->n_bg_frac, 0.2);
      RP::add<Real>(pop + "_Reconnection.nSpaceSamples",
              "Number of spatial sampling points per dimension", sP->nSpaceSamples, 1);
      RP::add<Real>(pop + "_Reconnection.nVelocitySamples",
              "Number of velocity sampling points per dimension", sP->nVelocitySamples, 1);
   }
}

// ---------------------------------------------------------------------------
void Reconnection::getParameters() {
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      this->speciesParams.push_back(*this->speciesParamsRead.at(i));
   }
}

// ---------------------------------------------------------------------------
// Number density profile: Harris sheet + uniform background
Real Reconnection::getDensity(creal z, const uint popID) const {
   const ReconnectionSpeciesParameters& sP = speciesParams[popID];
   const Real sech = 1.0 / std::cosh(z / lambda);
   return sP.n * (sech * sech + sP.n_bg_frac);
}

// ---------------------------------------------------------------------------
// fillPhaseSpace: called once per cell/block at initialization.
// Fills all velocity-space blocks that have been allocated for this cell.
Realf Reconnection::fillPhaseSpace(spatial_cell::SpatialCell* cell,
                                   const uint popID,
                                   const uint nRequested) const {
   const ReconnectionSpeciesParameters& sP = speciesParams[popID];
   const Real mass = getObjectWrapper().particleSpecies[popID].mass;

   const Real z       = cell->parameters[CellParams::ZCRD] + 0.5 * cell->parameters[CellParams::DZ];
   const Real initRho = getDensity(z, popID);
   const Real initT   = sP.T;
   // Bulk velocity is zero everywhere initially
   const Real V0X = 0.0, V0Y = 0.0, V0Z = 0.0;

#ifdef USE_GPU
   vmesh::VelocityMesh* vmesh = cell->dev_get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer* VBC = cell->dev_get_velocity_blocks(popID);
#else
   vmesh::VelocityMesh* vmesh = cell->get_velocity_mesh(popID);
   vmesh::VelocityBlockContainer* VBC = cell->get_velocity_blocks(popID);
#endif

   Realf rhosum = 0;
   arch::parallel_reduce<arch::null>(
      {WID, WID, WID, nRequested},
      ARCH_LOOP_LAMBDA(const uint i, const uint j, const uint k,
                       const uint initIndex, Realf* lsum) {
         vmesh::GlobalID* GIDlist = vmesh->getGrid()->data();
         Realf* bufferData        = VBC->getData();
         const vmesh::GlobalID blockGID = GIDlist[initIndex];

         Real blockCoords[6];
         vmesh->getBlockInfo(blockGID, &blockCoords[0]);
         creal vxBlock = blockCoords[0];
         creal vyBlock = blockCoords[1];
         creal vzBlock = blockCoords[2];
         creal dvxCell = blockCoords[3];
         creal dvyCell = blockCoords[4];
         creal dvzCell = blockCoords[5];

         ARCH_INNER_BODY(i, j, k, initIndex, lsum) {
            creal vx = vxBlock + (i + 0.5) * dvxCell - V0X;
            creal vy = vyBlock + (j + 0.5) * dvyCell - V0Y;
            creal vz = vzBlock + (k + 0.5) * dvzCell - V0Z;
            const Realf value = MaxwellianPhaseSpaceDensity(vx, vy, vz, initT, initRho, mass);
            bufferData[initIndex * WID3 + k * WID2 + j * WID + i] = value;
         };
      },
      rhosum);
   return rhosum;
}

// ---------------------------------------------------------------------------
// probePhaseSpace: used by TriAxisSearch to probe a single velocity point.
Realf Reconnection::probePhaseSpace(spatial_cell::SpatialCell* cell,
                                    const uint popID,
                                    Real vx_in, Real vy_in, Real vz_in) const {
   const ReconnectionSpeciesParameters& sP = speciesParams[popID];
   const Real mass = getObjectWrapper().particleSpecies[popID].mass;

   const Real z       = cell->parameters[CellParams::ZCRD] + 0.5 * cell->parameters[CellParams::DZ];
   const Real initRho = getDensity(z, popID);
   const Real initT   = sP.T;

   const Realf value = MaxwellianPhaseSpaceDensity(vx_in, vy_in, vz_in, initT, initRho, mass);
   return value;
}

// ---------------------------------------------------------------------------
void Reconnection::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

// ---------------------------------------------------------------------------
// getV0: bulk velocity hint for TriAxisSearch block allocation.
// Zero bulk velocity everywhere.
std::vector<std::array<Real, 3>> Reconnection::getV0(
   creal x, creal y, creal z, const uint popID) const {
   std::vector<std::array<Real, 3>> V0;
   V0.push_back({0.0, 0.0, 0.0});
   return V0;
}

// ---------------------------------------------------------------------------
// setProjectBField: set the initial magnetic field via the fsgrid stencil loop.
//
// Harris equilibrium (single current sheet at z = 0):
//   Bx(z) = Bx0 * tanh(z/lambda)
//   By    = By0  (guide field)
//   Bz    = 0
//
// GEM-challenge perturbation (flux function psi = Psi0*cos(2pi*x/Lx)*cos(pi*z/Lz)):
//   dBx = -dpsi/dz = +Psi0 * (pi/Lz)  * cos(2pi*x/Lx) * sin(pi*z/Lz)
//   dBz = +dpsi/dx = -Psi0 * (2pi/Lx) * sin(2pi*x/Lx) * cos(pi*z/Lz)
//
void Reconnection::setProjectBField(fsgrids::perbspan perb,
                                    fsgrids::bgbspan bgb,
                                    fsgrids::technicalspan technical, FieldSolverGrid& fsgrid) {
   setBackgroundFieldToZero(fsgrid, technical, bgb);

   if (P::isRestart) return;

   const Real Lx = P::xmax - P::xmin;
   const Real Lz = P::zmax - P::zmin;

   // local copies for lambda capture
   const auto Bx0_l   = this->Bx0;
   const auto By0_l   = this->By0;
   const auto Bz0_l   = this->Bz0;
   const auto Psi0_l  = this->Psi0;
   const auto lambda_l = this->lambda;

   fsgrid.parallel_for([](int timerId) -> phiprof::Timer { return phiprof::Timer{timerId}; },
                       phiprof::initializeTimer("setProjectBField-loop"), technical,
                       [=](const fsgrid::Coordinates& coordinates, const fsgrid::FsStencil& stencil,
                          cuint sysBoundaryFlag, cuint sysBoundaryLayer) {
      const std::array<Real, 3> xyz = coordinates.getPhysicalCoords(stencil.i, stencil.j, stencil.k);
      const std::array<Real, 3> gridSpacing = coordinates.physicalGridSpacing;
      auto& cell = perb[stencil.ooo()];

      // Cell-centre coordinates (face-centred field uses +0.5*d offsets)
      const Real xc = xyz[0] + 0.5 * gridSpacing[0];
      const Real zc = xyz[2] + 0.5 * gridSpacing[2];

      cell[fsgrids::bfield::PERBX] =
         Bx0_l * std::tanh(zc / lambda_l)
         + Psi0_l * (M_PI / Lz) * std::cos(2.0 * M_PI * xc / Lx) * std::sin(M_PI * zc / Lz);

      cell[fsgrids::bfield::PERBY] = By0_l;

      cell[fsgrids::bfield::PERBZ] =
         Bz0_l
         - Psi0_l * (2.0 * M_PI / Lx) * std::sin(2.0 * M_PI * xc / Lx) * std::cos(M_PI * zc / Lz);
   });
}

} // namespace projects


