/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Vlasiator. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "spatial_cell.hpp"
#include "common.h"
#include "project.h"
#include "parameters.h"
#include "readparameters.h"
#include "vlasovmover.h"

using namespace std;

typedef flowthroughParameters FTP;
Real FTP::T = NAN;
Real FTP::rho = NAN;
Real FTP::Bx = NAN;
Real FTP::By = NAN;
Real FTP::Bz = NAN;
uint FTP::nSpaceSamples = 0;
uint FTP::nVelocitySamples = 0;

bool initializeProject(void) {return true;}

bool addProjectParameters(){
   typedef Readparameters RP;
   RP::add("Flowthrough.rho", "Number density (m^-3)", 0.0);
   RP::add("Flowthrough.T", "Temperature (K)", 0.0);
   RP::add("Flowthrough.Bx", "Magnetic field x component (T)", 0.0);
   RP::add("Flowthrough.By", "Magnetic field y component (T)", 0.0);
   RP::add("Flowthrough.Bz", "Magnetic field z component (T)", 0.0);
   RP::add("Flowthrough.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Flowthrough.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   return true;
}

bool getProjectParameters(){
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
   typedef Readparameters RP;
   if(!RP::get("Flowthrough.rho", FTP::rho)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Flowthrough.T", FTP::T)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Flowthrough.Bx", FTP::Bx)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Flowthrough.By", FTP::By)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Flowthrough.Bz", FTP::Bz)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Flowthrough.nSpaceSamples", FTP::nSpaceSamples)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   if(!RP::get("Flowthrough.nVelocitySamples", FTP::nVelocitySamples)) {
      if(myRank == MASTER_RANK) cerr << __FILE__ << ":" << __LINE__ << " ERROR: This option has not been added!" << endl;
      exit(1);
   }
   return true;
}

void setProjectCell(SpatialCell* cell) {
   // Set up cell parameters:
   calcCellParameters(&((*cell).parameters[0]), 0.0);
   
   cell->parameters[CellParams::RHOLOSSADJUST] = 0.0;
   cell->parameters[CellParams::RHOLOSSVELBOUNDARY] = 0.0;
   
   // Go through each velocity block in the velocity phase space grid.
   // Set the initial state and block parameters:
   creal dvx_block = SpatialCell::block_dvx; // Size of a block in vx-direction
   creal dvy_block = SpatialCell::block_dvy; //                    vy
   creal dvz_block = SpatialCell::block_dvz; //                    vz
   creal dvx_blockCell = SpatialCell::cell_dvx; // Size of one cell in a block in vx-direction
   creal dvy_blockCell = SpatialCell::cell_dvy; //                                vy
   creal dvz_blockCell = SpatialCell::cell_dvz; //                                vz
   
   for (uint kv=0; kv<P::vzblocks_ini; ++kv) 
      for (uint jv=0; jv<P::vyblocks_ini; ++jv)
         for (uint iv=0; iv<P::vxblocks_ini; ++iv) {
            creal vx_block = P::vxmin + iv*dvx_block; // vx-coordinate of the lower left corner
            creal vy_block = P::vymin + jv*dvy_block; // vy-
            creal vz_block = P::vzmin + kv*dvz_block; // vz-
            
            // Calculate volume average of distrib. function for each cell in the block.
            for (uint kc=0; kc<WID; ++kc) 
               for (uint jc=0; jc<WID; ++jc) 
                  for (uint ic=0; ic<WID; ++ic) {
                     creal vx_cell = vx_block + ic*dvx_blockCell;
                     creal vy_cell = vy_block + jc*dvy_blockCell;
                     creal vz_cell = vz_block + kc*dvz_blockCell;
                     Real average = 
                     calcPhaseSpaceDensity(cell->parameters[CellParams::XCRD],
                                           cell->parameters[CellParams::YCRD],
                                           cell->parameters[CellParams::ZCRD],
                                           cell->parameters[CellParams::DX],
                                           cell->parameters[CellParams::DY],
                                           cell->parameters[CellParams::DZ],
                                           vx_cell,vy_cell,vz_cell,
                                           dvx_blockCell,dvy_blockCell,dvz_blockCell);
                     
                     if(average!=0.0){
                        creal vx_cell_center = vx_block + (ic+convert<Real>(0.5))*dvx_blockCell;
                        creal vy_cell_center = vy_block + (jc+convert<Real>(0.5))*dvy_blockCell;
                        creal vz_cell_center = vz_block + (kc+convert<Real>(0.5))*dvz_blockCell;
                        cell->set_value(vx_cell_center,vy_cell_center,vz_cell_center,average);
                     }
                  }
         }
   calculateCellVelocityMoments(cell);
   
   //let's get rid of blocks not fulfilling the criteria here to save memory.
   cell->adjustSingleCellVelocityBlocks();
}

Real getDistribValue(creal& x,creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   return FTP::rho * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * FTP::T), 1.5) *
   exp(- physicalconstants::MASS_PROTON * (vx*vx + vy*vy + vz*vz) / (2.0 * physicalconstants::K_B * FTP::T));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {   
   creal d_x = dx / (FTP::nSpaceSamples-1);
   creal d_y = dy / (FTP::nSpaceSamples-1);
   creal d_z = dz / (FTP::nSpaceSamples-1);
   creal d_vx = dvx / (FTP::nVelocitySamples-1);
   creal d_vy = dvy / (FTP::nVelocitySamples-1);
   creal d_vz = dvz / (FTP::nVelocitySamples-1);
   
   Real avg = 0.0;
// #pragma omp parallel for collapse(6) reduction(+:avg)
   // WARNING No threading here if calling functions are already threaded
   for (uint i=0; i<FTP::nSpaceSamples; ++i)
      for (uint j=0; j<FTP::nSpaceSamples; ++j)
         for (uint k=0; k<FTP::nSpaceSamples; ++k)
            for (uint vi=0; vi<FTP::nVelocitySamples; ++vi)
               for (uint vj=0; vj<FTP::nVelocitySamples; ++vj)
                  for (uint vk=0; vk<FTP::nVelocitySamples; ++vk) {
                     avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz, dvx, dvy, dvz);
                  }
   return avg / pow(FTP::nSpaceSamples, 3.0) / pow(FTP::nVelocitySamples, 3.0);
   
//    CellID cellID = 1 + round((x - Parameters::xmin) / dx + 
//    (y - Parameters::ymin) / dy * Parameters::xcells_ini +
//    (z - Parameters::zmin) / dz * Parameters::ycells_ini * Parameters::xcells_ini);
   
//    return cellID * pow(physicalconstants::MASS_PROTON / (2.0 * M_PI * physicalconstants::K_B * FTP::T), 1.5) *
//    exp(- physicalconstants::MASS_PROTON * (vx*vx + vy*vy + vz*vz) / (2.0 * physicalconstants::K_B * FTP::T));
}

bool cellParametersChanged(creal& t) {return false;}

void calcBlockParameters(Real* blockParams) { }

void calcCellParameters(Real* cellParams,creal& t) {
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = FTP::Bx;
   cellParams[CellParams::BY   ] = FTP::By;
   cellParams[CellParams::BZ   ] = FTP::Bz;
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->parameters, t);
   }
}

void dipole(creal x, creal y, creal z, Real& Bx, Real &By, Real& Bz) {
   creal k_0 = 8.0e15; // Wb m
   Real r = sqrt(x*x + y*y + z*z); // radial
   Real theta = atan2(sqrt(x*x + y*y), z); // polar
   Real phi = atan2(y, x); // azimuthal
   
   Bx = sin(theta) * cos(theta) * cos(phi);
   By = sin(theta) * cos(theta) * sin(phi);
   Bz = cos(theta)*cos(theta) - 1.0 / 3.0;
   
   Bx *= 3.0 * k_0 / (r*r*r);
   By *= 3.0 * k_0 / (r*r*r);
   Bz *= 3.0 * k_0 / (r*r*r);
}
