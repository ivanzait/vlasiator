#ifndef CPU_ACC_H
#define CPU_ACC_H

#include "definitions.h"
#include "common.h"
#include "cell_spatial.h"
#include "project.h"
#include "leveque_common.h"

// Constant for switch statement in solver to decide which of the eight 
// possibilities should be calculated. 
cuint AXP_AYP_AZP = 0; // Ax > 0, Ay > 0, Az > 0
cuint AXP_AYP_AZN = 1; //                 Az < 0
cuint AXP_AYN_AZP = 2; //         Ay < 0, Az > 0
cuint AXP_AYN_AZN = 3; //                 Az < 0
cuint AXN_AYP_AZP = 4; // Ax < 0, Ay > 0, Az > 0
cuint AXN_AYP_AZN = 5; //                 Az < 0
cuint AXN_AYN_AZP = 6; //         Ay < 0, Az > 0
cuint AXN_AYN_AZN = 7; //                 Az < 0

template<typename T> inline T accIndex(const T& i,const T& j,const T& k) {return k*WID2+j*WID+i;}
template<typename T> inline T fullInd(const T& i,const T& j,const T& k) {return k*64+j*8+i;}
template<typename T> inline T isBoundary(const T& STATE,const T& BND) {return STATE & BND;}
template<typename T> inline T findex(const T& i,const T& j,const T& k) {return k*36+j*6+i;}

template<typename REAL,typename UINT> void accumulateChanges(const UINT& BLOCK,const REAL* const dF,REAL* const flux,const UINT* const nbrsVel) {
   UINT nbrBlock;
   REAL* nbrFlux;
   const UINT STATE = nbrsVel[NbrsVel::STATE];
   typedef Parameters P;

   // NOTE: velocity block can have up to 26 neighbours, and we need to copy changes 
   // to each existing neighbour here.
   // NOTE2: dF is a (6,6,6) sized block, nbrFlux is (4,4,4).
   // NOTE3: This function copies values correctly
   const UINT MIN = 0;
   const UINT MAX = 5;
   const UINT ACCMIN = 0;
   const UINT ACCMAX = 3;
   
   // Accumulate changes to this block:
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) 
     flux[BLOCK*SIZE_FLUXS+accIndex(i,j,k)] += dF[findex(i+1,j+1,k+1)];
   
   // Accumulate changes to (iv-1,jv,kv) neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VXNEG];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,k)] += dF[findex(MIN,j+1,k+1)];
      // Accumulate changes to (iv-1,jv+-1,kv) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMAX,ACCMAX,k)] += dF[findex(MIN,MIN,k+1)];
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMAX,ACCMIN,k)] += dF[findex(MIN,MAX,k+1)];
      }
   }
   // Accumulate changes to (iv+1,jv,kv) neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VXPOS];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,k)] += dF[findex(MAX,j+1,k+1)];
      // Accumulate changes to (iv+1,jv+-1,kv) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMIN,ACCMAX,k)] += dF[findex(MAX,MIN,k+1)];
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT k=0; k<WID; ++k) nbrFlux[accIndex(ACCMIN,ACCMIN,k)] += dF[findex(MAX,MAX,k+1)];
      }
   }
   // Accumulate changes to -vy neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VYNEG];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,k)] += dF[findex(i+1,MIN,k+1)];
      // Accumulate changes to (iv,jv-1,kv+-1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,ACCMAX)] += dF[findex(i+1,MIN,MIN)];
      }
      if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMAX,ACCMIN)] += dF[findex(i+1,MIN,MAX)];
      }
   }
   // Accumulate changes to +vy neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VYPOS];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT k=0; k<WID; ++k) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,k)] += dF[findex(i+1,MAX,k+1)];
      // Accumulate changes to (iv,jv+1,kv+-1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini; // temp solution
      	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,ACCMAX)] += dF[findex(i+1,MAX,MIN)];
      }
      if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,ACCMIN,ACCMIN)] += dF[findex(i+1,MAX,MAX)];
      }
   }
   // Accumulate changes to -vz neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VZNEG];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,j,ACCMAX)] += dF[findex(i+1,j+1,MIN)];
      // Accumulate changes to (iv+-1,jv,kv-1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,ACCMAX)] += dF[findex(MIN,j+1,MIN)];
      }
      if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
	 nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,ACCMAX)] += dF[findex(MAX,j+1,MIN)];
      }
   }
   // Accumulate changes to +vz neighbour if it exists:
   if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) {
      nbrBlock = nbrsVel[NbrsVel::VZPOS];
      nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
      for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) nbrFlux[accIndex(i,j,ACCMIN)] += dF[findex(i+1,j+1,MAX)];
      // Accumulate changes to (iv+-1,jv,kv+1) neighbours if they exist:
      if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMAX,j,ACCMIN)] += dF[findex(MIN,j+1,MAX)];
      }
      if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
	 nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + 1; // temp solution
	 nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	 for (UINT j=0; j<WID; ++j) nbrFlux[accIndex(ACCMIN,j,ACCMIN)] += dF[findex(MAX,j+1,MAX)];
      }
   }
   
   // Accumulate changes to 8 corner neighbours:
   if (isBoundary(STATE,NbrsVel::VX_NEG_BND) == 0) {
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv-1,jv-1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMAX,ACCMAX)] += dF[findex(MIN,MIN,MIN)];
	 }
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv-1,jv-1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMAX,ACCMIN)] += dF[findex(MIN,MIN,MAX)];
	 }
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) {
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv-1,jv+1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMIN,ACCMAX)] += dF[findex(MIN,MAX,MIN)];
	 } 
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv-1,jv+1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini - 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMAX,ACCMIN,ACCMIN)] += dF[findex(MIN,MAX,MAX)];
	 }
      }
   }
   if (isBoundary(STATE,NbrsVel::VX_POS_BND) == 0) {
      if (isBoundary(STATE,NbrsVel::VY_NEG_BND) == 0) {
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv+1,jv-1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMAX,ACCMAX)] += dF[findex(MAX,MIN,MIN)];
	 }
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv+1,jv-1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini - P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMAX,ACCMIN)] += dF[findex(MAX,MIN,MAX)];
	 }
      }
      if (isBoundary(STATE,NbrsVel::VY_POS_BND) == 0) { 
	 if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) == 0) { // (iv+1,jv+1,kv-1)
	    nbrBlock = BLOCK - P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMIN,ACCMAX)] += dF[findex(MAX,MAX,MIN)];
	 }
	 if (isBoundary(STATE,NbrsVel::VZ_POS_BND) == 0) { // (iv+1,jv+1,kv+1)
	    nbrBlock = BLOCK + P::vyblocks_ini*P::vxblocks_ini + P::vxblocks_ini + 1; // temp solution
	    nbrFlux  = flux + nbrBlock*SIZE_FLUXS;
	    nbrFlux[accIndex(ACCMIN,ACCMIN,ACCMIN)] += dF[findex(MAX,MAX,MAX)];
	 }
      }
   }
}

template<typename REAL,typename UINT> void fetchAllAverages(const UINT& BLOCK,REAL* const avgs,const REAL* const cpu_avgs,const UINT* const nbrsVel) {
   UINT nbrBlock;
   for (UINT i=0; i<8*WID3; ++i) avgs[i] = 0.0;
   const UINT STATE = nbrsVel[BLOCK*SIZE_NBRS_VEL+NbrsVel::STATE];

   // Copy averages from -x neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VX_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<2; ++i) {
	 avgs[fullInd(i  ,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK; 
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<2; ++i) {
	 //avgs[fullInd(i  ,j+2,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i+2,j,k)];
	 avgs[fullInd(i  ,j+2,k+2)] = tmp[accIndex(i+2,j,k)];
	 }
      }
   }
   // Copy averages from +x neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VX_POS_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<2; ++i) {
	 avgs[fullInd(i+6,j+2,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VXPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<2; ++i) {
	    //avgs[fullInd(i+6,j+2,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i  ,j,k)];
	    avgs[fullInd(i+6,j+2,k+2)] = tmp[accIndex(i  ,j,k)];
	 }
      }
   }
   // Copy averages from -y neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VY_NEG_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j  ,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j  ,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j+2,k)];
	    avgs[fullInd(i+2,j  ,k+2)] = tmp[accIndex(i,j+2,k)];
	 }
      }
   }
   // Copy averages from +y neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VY_POS_BND) > 0) {
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+6,k+2)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VYPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<WID; ++k) for (UINT j=0; j<2; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j+6,k+2)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j  ,k)];
	    avgs[fullInd(i+2,j+6,k+2)] = tmp[accIndex(i,j  ,k)];
	 }
      }
   }
   // Copy averages from -z neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VZ_NEG_BND) > 0) {
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k  )] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZNEG];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j+2,k  )] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,k+2)];
	    avgs[fullInd(i+2,j+2,k  )] = tmp[accIndex(i,j,k+2)];
	 }
      }
   }
   // Copy averages from +z neighbour, or calculate using a boundary function:
   if (isBoundary(STATE,NbrsVel::VZ_POS_BND) > 0) {
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
	 avgs[fullInd(i+2,j+2,k+6)] = 0.0; // BOUNDARY VALUE
      }
   } else {
      nbrBlock = nbrsVel[BLOCK*SIZE_NBRS_VEL + NbrsVel::VZPOS];
      creal* const tmp = cpu_avgs + nbrBlock*SIZE_VELBLOCK;
      for (UINT k=0; k<2; ++k) for (UINT j=0; j<WID; ++j) {
	 #pragma ivdep
	 for (UINT i=0; i<WID; ++i) {
	    //avgs[fullInd(i+2,j+2,k+6)] = cpu_avgs[nbrBlock*SIZE_VELBLOCK + accIndex(i,j,k)];
	    avgs[fullInd(i+2,j+2,k+6)] = tmp[accIndex(i,j,k)];
	 }
      }
   }

   // Copy volume averages of this block:
   creal* const tmp = cpu_avgs + BLOCK*SIZE_VELBLOCK;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      //avgs[fullInd(i+2,j+2,k+2)] = cpu_avgs[BLOCK*SIZE_VELBLOCK + accIndex(i,j,k)];
      avgs[fullInd(i+2,j+2,k+2)] = tmp[accIndex(i,j,k)];
   }
}

template<typename REAL,typename UINT,typename CELL> void cpu_clearVelFluxes(CELL& cell,const UINT& BLOCK) {
   for (UINT i=0; i<SIZE_FLUXS; ++i) cell.cpu_fx[BLOCK*SIZE_FLUXS + i] = 0.0;
}

template<typename REAL,typename UINT,typename CELL> void cpu_calcVelFluxes(CELL& cell,const UINT& BLOCK,const REAL& DT,creal* const accmat) {
   // Creation of temporary calculation block dfdt and avgs + 
   // value fetching and initializations seem to take about
   // ~4% of time used by calcVelFluxes
   
   // Allocate temporary array in which local dF changes are calculated:
   // F is the distribution function, dFdt is its change over timestep DT
   Real dfdt[216];
   for (UINT i=0; i<216; ++i) dfdt[i] = 0.0;

   const REAL* const cellParams = cell.cpu_cellParams;
   const REAL* const blockParams = cell.cpu_blockParams + BLOCK*SIZE_BLOCKPARAMS; 
   REAL avgs[8*WID3];
   fetchAllAverages(BLOCK,avgs,cell.cpu_avgs,cell.cpu_nbrsVel);

   const REAL DVX = blockParams[BlockParams::DVX];
   const REAL DVY = blockParams[BlockParams::DVY];
   const REAL DVZ = blockParams[BlockParams::DVZ];
   const REAL dt_per_dvx = DT / blockParams[BlockParams::DVX];
   const REAL dt_per_dvy = DT / blockParams[BlockParams::DVY];
   const REAL dt_per_dvz = DT / blockParams[BlockParams::DVZ];
   
   REAL Ax,Ay,Az;
   for (UINT k=0; k<WID; ++k) for (UINT j=0; j<WID; ++j) for (UINT i=0; i<WID; ++i) {
      REAL R,incrWave,transIncrWave,doubleTransIncrWave;
      REAL corrWave,transCorrWave,doubleTransCorrWave;
      UINT solverFlags;
      
      // ***********************************
      // ***** INTERFACE BETWEEN I,I-1 *****
      // ***********************************
      const REAL xcc = avgs[fullInd(i+2,j+2,k+2)];
      const REAL xp1 = avgs[fullInd(i+3,j+2,k+2)];
      const REAL xm1 = avgs[fullInd(i+1,j+2,k+2)];
      const REAL xm2 = avgs[fullInd(i  ,j+2,k+2)];
      calcAccFaceX(Ax,Ay,Az,i,j,k,cellParams,blockParams);

      solverFlags = 0;
      if (Az < ZERO) solverFlags = (solverFlags | (1 << 0));
      if (Ay < ZERO) solverFlags = (solverFlags | (1 << 1));
      if (Ax < ZERO) solverFlags = (solverFlags | (1 << 2));

      switch (solverFlags) {
       case AXP_AYP_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
	 #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;

	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+2)] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvx*Ax)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYP_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;

	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k  )] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 #endif
	 break;
       case AXP_AYN_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;

	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+2)] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE-Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xm1*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ax*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
         #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvx*Ax)*R*limiter(xm1-xm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE-Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xm1-xm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
	 #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+2)] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+2,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+2)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] += doubleTransCorrWave;
         #endif
	 break;	 
       case AXN_AYP_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
         dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k  )] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
         #endif
	 break;	 
       case AXN_AYN_AZP:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += transCorrWave;
	 #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+2)] += TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZN:
	 R = xcc-xm1;
	 incrWave = Ax*xcc*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ax*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvx*Ax)*R*limiter(xp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
         #endif

	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k  )] -= TWO*doubleTransIncrWave*R; // x
	 dfdt[findex(i  ,j+1,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE+Ax*dt_per_dvx)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(xp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
      }
      
      // ***********************************
      // ***** INTERFACE BETWEEN J,J-1 *****
      // ***********************************
      const REAL yp1 = avgs[fullInd(i+2,j+3,k+2)];
      const REAL ym1 = avgs[fullInd(i+2,j+1,k+2)];
      const REAL ym2 = avgs[fullInd(i+2,j  ,k+2)];
      calcAccFaceY(Ax,Ay,Az,i,j,k,cellParams,blockParams);
      
      solverFlags = 0;
      if (Az < ZERO) solverFlags = (solverFlags | (1 << 0));
      if (Ay < ZERO) solverFlags = (solverFlags | (1 << 1));
      if (Ax < ZERO) solverFlags = (solverFlags | (1 << 2));

      switch (solverFlags) {
       case AXP_AYP_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
	 #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+2)] -= TWO*doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] += doubleTransCorrWave;
	 #endif
	 break;
       case AXP_AYP_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k  )] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
	 #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+2)] -= TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j  ,k+2)] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j  ,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] += doubleTransCorrWave;
	 #endif	 
	 break;
       case AXP_AYN_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k  )] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
	 #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+2)] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j+1,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*ym1*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Ay*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvy*Ay)*R*limiter(ym1-ym2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] -= TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j+1,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(ym1-ym2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZP:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+2)] -= transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+2)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+2)] += TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+2)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+2)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+2)] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZN:
	 R = xcc-ym1;
	 incrWave = Ay*xcc*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= corrWave;
         #endif

	 transIncrWave = HALF*Ax*Ay*dt_per_dvx*dt_per_dvy - SIXTH*Ax*Ay*fabs(Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j  ,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Ay*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvy;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvy*Ay)*R*limiter(yp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k  )] -= TWO*doubleTransIncrWave*R; // y
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvy*Ay)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(yp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
         #endif
	 break;
      }

      // ***********************************
      // ***** INTERFACE BETWEEN K,K-1 *****
      // ***********************************
      const REAL zp1 = avgs[fullInd(i+2,j+2,k+3)];
      const REAL zm1 = avgs[fullInd(i+2,j+2,k+1)];
      const REAL zm2 = avgs[fullInd(i+2,j+2,k  )];
      calcAccFaceZ(Ax,Ay,Az,i,j,k,cellParams,blockParams);
      
      solverFlags = 0;
      if (Az < ZERO) solverFlags = (solverFlags | (1 << 0));
      if (Ay < ZERO) solverFlags = (solverFlags | (1 << 1));
      if (Ax < ZERO) solverFlags = (solverFlags | (1 << 2));
      
      switch (solverFlags) {
       case AXP_AYP_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
	 #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+2,k+1)] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] += doubleTransCorrWave;
	 #endif
	 break;
       case AXP_AYP_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+2,k  )] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j+2,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k+1)] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXP_AYN_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+2,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+2,j  ,k  )] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+2,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+2,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+2,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k+1)] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransIncrWave*R;
	 dfdt[findex(i+1,j+2,k+1)] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYP_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k  )] -= transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] += transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+2,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+2,k  )] += TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+2,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+2,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZP:
	 R = xcc - zm1;
	 incrWave = Az*zm1*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
         #if MOVER_VLASOV_ORDER > 1
	 corrWave = HALF*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k+1)] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k+1)] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = HALF*Ax*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = HALF*Ay*Az*(ONE - dt_per_dvz*Az)*R*limiter(zm1-zm2,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k+1)] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i+1,j  ,k+1)] += doubleTransIncrWave*R;
	 dfdt[findex(i  ,j+1,k+1)] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = HALF*Ax*Ay*Az*(ONE - dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zm1-zm2,R+EPSILON,xcc);
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
       case AXN_AYN_AZN:
	 R = xcc - zm1;
	 incrWave = Az*xcc*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += incrWave;
	 dfdt[findex(i+1,j+1,k  )] -= incrWave;
	 #if MOVER_VLASOV_ORDER > 1
	 corrWave = -HALF*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += corrWave;
	 dfdt[findex(i+1,j+1,k  )] -= corrWave;
         #endif
	 
	 transIncrWave = HALF*Ax*Az*dt_per_dvx*dt_per_dvz - SIXTH*Ax*fabs(Ay)*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j+1,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 transIncrWave = HALF*Ay*Az*dt_per_dvy*dt_per_dvz - SIXTH*fabs(Ax)*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j  ,k  )] += transIncrWave*R;
	 dfdt[findex(i+1,j+1,k  )] -= transIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 transCorrWave = -HALF*Ax*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvx*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += transCorrWave;
	 transCorrWave = -HALF*Ay*Az*(ONE + dt_per_dvz*Az)*R*limiter(zp1-xcc,R+EPSILON,xcc)*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i+1,j+1,k+1)] += transCorrWave;
	 dfdt[findex(i+1,j  ,k+1)] -= transCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= transCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += transCorrWave;
         #endif
	 
	 doubleTransIncrWave = SIXTH*Ax*Ay*Az*dt_per_dvx*dt_per_dvy*dt_per_dvz;
	 dfdt[findex(i  ,j  ,k  )] -= TWO*doubleTransIncrWave*R; // z
	 dfdt[findex(i  ,j+1,k  )] += doubleTransIncrWave*R;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransIncrWave*R;
	 #if MOVER_VLASOV_ORDER > 1
	 doubleTransCorrWave = -HALF*Ax*Ay*Az*(ONE + dt_per_dvz*Az)*dt_per_dvx*dt_per_dvy*dt_per_dvz*R*limiter(zp1-xcc,R+EPSILON,xcc);
	 dfdt[findex(i+1,j  ,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i+1,j  ,k  )] += doubleTransCorrWave;
	 dfdt[findex(i  ,j  ,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k+1)] += doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k+1)] -= doubleTransCorrWave;
	 dfdt[findex(i+1,j+1,k  )] -= doubleTransCorrWave;
	 dfdt[findex(i  ,j+1,k  )] += doubleTransCorrWave;
         #endif
	 break;
      }

   }
   
   // If multithreading is used (OpenMP/pthreads), then some of the updates 
   // in accumulateChanges may need to be atomic. dfdt values are copied to 
   // neighbouring blocks within the same spatial cell.
   // 
   // accumulateChanges seems to take ~3% of time used by calcVelFluxes
   accumulateChanges(BLOCK,dfdt,cell.cpu_fx,cell.cpu_nbrsVel+BLOCK*SIZE_NBRS_VEL);
}

template<typename REAL,typename UINT,typename CELL> void cpu_propagateVel(CELL& cell,const UINT& BLOCK,const REAL& DT) {
   for (UINT i=0; i<WID3; ++i) cell.cpu_avgs[BLOCK*SIZE_VELBLOCK+i] += cell.cpu_fx[BLOCK*SIZE_VELBLOCK+i];
}
#endif
