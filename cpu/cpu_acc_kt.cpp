#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "cpu_acc_kt.h"

using namespace std;

bool cpu_acceleration(SpatialCell& cell) {
   // Clear spatial cell velocity moments:
   cell.cpu_cellParams[CellParams::RHO]   = 0.0;
   cell.cpu_cellParams[CellParams::RHOVX] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVY] = 0.0;
   cell.cpu_cellParams[CellParams::RHOVZ] = 0.0;
   
   bool success = true;
   const Real DT = Parameters::dt;
   
   // First pass of vx-propagation
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelX<Real>(cell,block,0.5*DT);
   
   // First pass of vy-propagation
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsY<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesY<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelY<Real>(cell,block,DT);
   /*
   // vz-propagation
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsZ<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesZ<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelZ<Real>(cell,block,DT);
   
   // Second pass of vy-propagation
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsY<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesY<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelY<Real>(cell,block,0.5*DT);
   */
   // Second pass of vx-propagation
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelDerivsX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_calcVelFluxesX<Real>(cell,block);
   for (uint block=0; block<cell.N_blocks; ++block) cpu_propagateVelX<Real>(cell,block,0.5*DT);
   
   return success;
}

