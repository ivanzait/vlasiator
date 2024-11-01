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

#include <cstdlib>
#include <iostream>
#include <iomanip> // for setprecision()
#include <cmath>
#include <vector>
#include <sstream>
#include <ctime>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "grid.h"
#include "vlasovsolver/vlasovmover.h"
#include "definitions.h"
#include "mpiconversion.h"
#include "logger.h"
#include "parameters.h"
#include "datareduction/datareducer.h"
#include "sysboundary/sysboundary.h"
#include "fieldsolver/fs_common.h"
#include "fieldsolver/gridGlue.hpp"
#include "fieldsolver/derivatives.hpp"
#include "vlasovsolver/cpu_trans_pencils.hpp"
#include "projects/project.h"
#include "iowrite.h"
#include "ioread.h"
#include "object_wrapper.h"

#ifdef PAPI_MEM
#include "papi.h"
#endif

#include "vlasovsolver/cpu_trans_pencils.hpp"
#ifdef USE_GPU
#include "arch/gpu_base.hpp"
#endif

#ifdef DEBUG_VLASIATOR
#define DEBUG_GRID
#endif

using namespace std;

int globalflags::AMRstencilWidth = VLASOV_STENCIL_WIDTH;

extern Logger logFile, diagnostic;

void initVelocityGridGeometry(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
void initSpatialCellCoordinates(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);
void initializeStencils(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

void writeVelMesh(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   const vector<CellID>& cells = getLocalCells();

   static int counter=0;

      stringstream fname;
   fname << "VelMesh.";
   fname.width(3);
   fname.fill(0);
   fname << counter << ".vlsv";

   vlsv::Writer vlsvWriter;
   vlsvWriter.open(fname.str(),MPI_COMM_WORLD,0,MPI_INFO_NULL);
   writeVelocityDistributionData(vlsvWriter,mpiGrid,cells,MPI_COMM_WORLD);
   vlsvWriter.close();

   ++counter;
}

void initializeGrids(
   int argn,
   char **argc,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   FsGrid< std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> & perBGrid,
   FsGrid< std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH> & BgBGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsGrid,
   FsGrid< std::array<Real, fsgrids::moments::N_MOMENTS>, FS_STENCIL_WIDTH> & momentsDt2Grid,
   FsGrid< std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> & EGrid,
   FsGrid< std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> & EGradPeGrid,
   FsGrid< std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> & volGrid,
   FsGrid< fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid,
   SysBoundary& sysBoundaries,
   Project& project
) {
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

   // Init Zoltan:
   float zoltanVersion;
   if (Zoltan_Initialize(argn,argc,&zoltanVersion) != ZOLTAN_OK) {
      if(myRank == MASTER_RANK) cerr << "\t ERROR: Zoltan initialization failed." << endl;
      exit(1);
   } else {
      logFile << "\t Zoltan " << zoltanVersion << " initialized successfully" << std::endl << writeVerbose;
   }

   MPI_Comm comm = MPI_COMM_WORLD;
   int neighborhood_size = VLASOV_STENCIL_WIDTH;
   if (P::amrMaxSpatialRefLevel > 0) {
      switch (VLASOV_STENCIL_WIDTH) {
         case 1:
            // Required cells will be included already
            break;
         case 2:
            // looking from high to low refinement: stencil 2 will only give 1 cell, so need to add 1
            neighborhood_size = VLASOV_STENCIL_WIDTH+1;
            break;
         case 3:
            // looking from high to low refinement: stencil 3 will only give 2 cells, so need to add 2
            // to reach surely into the third low-refinement neighbour
            neighborhood_size = VLASOV_STENCIL_WIDTH+2;
            break;
         default:
            std::cerr<<"Warning: unrecognized VLASOV_STENCIL_WIDTH in grid.cpp"<<std::endl;
      }
   }
   globalflags::AMRstencilWidth = neighborhood_size;

   const std::array<uint64_t, 3> grid_length = {{P::xcells_ini, P::ycells_ini, P::zcells_ini}};
   dccrg::Cartesian_Geometry::Parameters geom_params;
   geom_params.start[0] = P::xmin;
   geom_params.start[1] = P::ymin;
   geom_params.start[2] = P::zmin;
   geom_params.level_0_cell_length[0] = P::dx_ini;
   geom_params.level_0_cell_length[1] = P::dy_ini;
   geom_params.level_0_cell_length[2] = P::dz_ini;

   phiprof::Timer dccrgTimer {"Initialize DCCRG grid"};
   mpiGrid.set_initial_length(grid_length)
      .set_load_balancing_method(&P::loadBalanceAlgorithm[0])
      .set_neighborhood_length(neighborhood_size)
      .set_maximum_refinement_level(P::amrMaxSpatialRefLevel)
      .set_periodic(sysBoundaries.isPeriodic(0),
                    sysBoundaries.isPeriodic(1),
                    sysBoundaries.isPeriodic(2))
      .initialize(comm)
      .set_geometry(geom_params);
   dccrgTimer.stop();

   phiprof::Timer refineTimer {"Refine spatial cells"};
   // We need this first as well
   recalculateLocalCellsCache();
   if (!P::isRestart) {
      if (P::amrMaxSpatialRefLevel > 0 && project.refineSpatialCells(mpiGrid)) {
         mpiGrid.balance_load();
         recalculateLocalCellsCache();
         mapRefinement(mpiGrid, technicalGrid);
      }
   } else {
      if (readFileCells(mpiGrid, P::restartFileName)) {
         mpiGrid.balance_load();
         recalculateLocalCellsCache();
         mapRefinement(mpiGrid, technicalGrid);
      }
   }
   refineTimer.stop();
   initializeStencils(mpiGrid);

   for (const auto& [key, value] : P::loadBalanceOptions) {
      mpiGrid.set_partitioning_option(key, value);
   }
   phiprof::Timer initialLBTimer {"Initial load-balancing"};
   if (myRank == MASTER_RANK) logFile << "(INIT): Starting initial load balance." << endl << writeVerbose;
   mpiGrid.balance_load(); // Direct DCCRG call, recalculate cache afterwards
   recalculateLocalCellsCache();

   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
   mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);

   setFaceNeighborRanks( mpiGrid );
   const vector<CellID>& cells = getLocalCells();
   initialLBTimer.stop();

   if (myRank == MASTER_RANK) {
      logFile << "(INIT): Set initial state." << endl << writeVerbose;
   }

   phiprof::Timer initialStateTimer {"Set initial state"};

   phiprof::Timer setCoordsTimer {"Set spatial cell coordinates"};
   initSpatialCellCoordinates(mpiGrid);
   setCoordsTimer.stop();

   sysBoundaries.initSysBoundaries(project, P::t_min);

   SpatialCell::set_mpi_transfer_type(Transfer::CELL_DIMENSIONS);
   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_NEIGHBORHOOD_ID);

   // We want this before restart refinement
   phiprof::Timer classifyTimer {"Classify cells (sys boundary conditions)"};
   sysBoundaries.classifyCells(mpiGrid,technicalGrid);
   classifyTimer.stop();

   if (P::isRestart) {
      logFile << "Restart from "<< P::restartFileName << std::endl << writeVerbose;
      phiprof::Timer restartReadTimer {"Read restart"};
      if (readGrid(mpiGrid,perBGrid,EGrid,technicalGrid,P::restartFileName) == false) {
         logFile << "(MAIN) ERROR: restarting failed" << endl;
         exit(1);
      }
      restartReadTimer.stop();

      if (P::forceRefinement) {
         // Adapt refinement to match new static refinement parameters
         phiprof::Timer timer {"Restart refinement"};
         for (int i = 0; i < P::amrMaxSpatialRefLevel; ++i) {
            // (un)Refinement is done one level at a time so we don't blow up memory
            if (!adaptRefinement(mpiGrid, technicalGrid, sysBoundaries, project, i)) {
               cerr << "(MAIN) ERROR: Forcing refinement takes too much memory" << endl;
               exit(1);
            }
            balanceLoad(mpiGrid, sysBoundaries);
         }
      } else if (P::refineOnRestart) {
         // Considered deprecated
         phiprof::Timer timer {"Restart refinement"};
         // Get good load balancing for refinement
         balanceLoad(mpiGrid, sysBoundaries);
         adaptRefinement(mpiGrid, technicalGrid, sysBoundaries, project);
         balanceLoad(mpiGrid, sysBoundaries);
      }
   }

   // Check refined cells do not touch boundary cells
   phiprof::Timer boundaryCheckTimer {"Check boundary refinement"};
   sysBoundaries.checkRefinement(mpiGrid);
   boundaryCheckTimer.stop();

   #ifdef USE_GPU
   phiprof::Timer prefetchDeviceTimer {"prefetch to GPU"};
   for (size_t i=0; i<cells.size(); ++i) {
      SpatialCell* cell = mpiGrid[cells[i]];
      cell->prefetchDevice();
   }
   prefetchDeviceTimer.stop();
   #endif

   if (P::isRestart) {
      //initial state for sys-boundary cells, will skip those not set to be reapplied at restart
      sysBoundaries.applyInitialState(mpiGrid, technicalGrid, perBGrid, BgBGrid, project);
   }

   // Update technicalGrid
   technicalGrid.updateGhostCells(); // This needs to be done at some point

   if (!P::isRestart && !P::writeFullBGB) {
      // If we are starting a new regular simulation, we need to prepare all cells with their initial state.
      // If we're only after writing out the full BGB we don't need all this shebang EXCEPT the weights!

      //Initial state based on project, background field in all cells
      //and other initial values in non-sysboundary cells
      phiprof::Timer applyInitialTimer {"Apply initial state"};
      // Go through every cell on this node and initialize the
      //  -Background field on all cells
      //  -Perturbed fields and ion distribution function in non-sysboundary cells
      // Each initialization has to be independent to avoid threading problems

      // Allow the project to set up data structures for it's setCell calls
      project.setupBeforeSetCell(cells);

      phiprof::Timer setCellTimer {"setCell"};
      #pragma omp parallel for schedule(dynamic)
      for (size_t i=0; i<cells.size(); ++i) {
         SpatialCell* cell = mpiGrid[cells[i]];
         if (cell->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            project.setCell(cell);
         }
      }
      setCellTimer.stop();

      // Initial state for sys-boundary cells
      sysBoundaries.applyInitialState(mpiGrid, technicalGrid, perBGrid, BgBGrid, project);

      #pragma omp parallel for schedule(static)
      for (size_t i=0; i<cells.size(); ++i) {
         mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER] = 0;
      }

      for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         adjustVelocityBlocks(mpiGrid,cells,true,popID);
         // set initial LB metric based on number of blocks, all others
         // will be based on time spent in acceleration
         #pragma omp parallel for schedule(static)
         for (size_t i=0; i<cells.size(); ++i) {
            mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER] += mpiGrid[cells[i]]->get_number_of_velocity_blocks(popID);
         }
      }

      shrink_to_fit_grid_data(mpiGrid); //get rid of excess data already here

      /*
      // Apply boundary conditions so that we get correct initial moments
      sysBoundaries.applySysBoundaryVlasovConditions(mpiGrid,Parameters::t);

      //compute moments, and set them  in RHO* and RHO_*_DT2. If restart, they are already read in
      phiprof::Timer initMomentsTimer {"Init moments"};
      calculateInitialVelocityMoments(mpiGrid);
      initMomentsTimer.stop();
      */

   } else if (P::writeFullBGB) {
      // If, instead of starting a regular simulation, we are only writing out the background field, it is enough to set a dummy load balance value of 1 here.
      for (size_t i=0; i<cells.size(); ++i) {
         mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER] = 1;
      }
   }


   // Balance load before we transfer all data below
   balanceLoad(mpiGrid, sysBoundaries);
   // Function includes re-calculation of local cells cache

   phiprof::Timer fetchNeighbourTimer {"Fetch Neighbour data", {"MPI"}};
   // update complete cell spatial data for full stencil (
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_copies_of_remote_neighbors(FULL_NEIGHBORHOOD_ID);
   fetchNeighbourTimer.stop();

   phiprof::Timer setBTimer {"project.setProjectBField"};
   project.setProjectBField(perBGrid, BgBGrid, technicalGrid);
   setBTimer.stop();
   phiprof::Timer fsGridGhostTimer {"fsgrid-ghost-updates"};
   perBGrid.updateGhostCells();
   BgBGrid.updateGhostCells();
   EGrid.updateGhostCells();

   // This will only have the BGB set up properly at this stage but we need the BGBvol for the Vlasov boundaries below.
   volGrid.updateGhostCells();
   fsGridGhostTimer.stop();
   phiprof::Timer getFieldsTimer {"getFieldsFromFsGrid"};
   getFieldsFromFsGrid(volGrid, BgBGrid, EGradPeGrid, technicalGrid, mpiGrid, cells);
   getFieldsTimer.stop();

   setBTimer.stop();

   // If we only want the full BGB for writeout, we have it now and we can return early.
   if(P::writeFullBGB == true) {
      return;
   }

   if (P::isRestart == false) {
      // Apply boundary conditions so that we get correct initial moments
      sysBoundaries.applySysBoundaryVlasovConditions(mpiGrid,Parameters::t, true); // It doesn't matter here whether we put _R or _V moments

      //compute moments, and set them  in RHO* and RHO_*_DT2. If restart, they are already read in
      phiprof::Timer timer {"Init moments"};
      calculateInitialVelocityMoments(mpiGrid);
   } else {
      phiprof::Timer timer {"Init moments"};
      for (size_t i=0; i<cells.size(); ++i) {
         calculateCellMoments(mpiGrid[cells[i]], true, true);
      }
   }

   phiprof::Timer finishFSGridTimer {"Finish fsgrid setup"};
   feedMomentsIntoFsGrid(mpiGrid, cells, momentsGrid,technicalGrid, false);
   if(!P::isRestart) {
      // WARNING this means moments and dt2 moments are the same here at t=0, which is a feature so far.
      feedMomentsIntoFsGrid(mpiGrid, cells, momentsDt2Grid, technicalGrid, false);
   } else {
      feedMomentsIntoFsGrid(mpiGrid, cells, momentsDt2Grid, technicalGrid, true);
   }
   momentsGrid.updateGhostCells();
   momentsDt2Grid.updateGhostCells();
   finishFSGridTimer.stop();

   // Set this so CFL doesn't break
   if(P::refineOnRestart) {
      // Half-step acceleration
      if( P::propagateVlasovAcceleration ) {
         calculateAcceleration(mpiGrid, -0.5*P::dt + 0.5*P::bailout_min_dt);
      } else {
         calculateAcceleration(mpiGrid, 0.0);
      }
      P::dt = P::bailout_min_dt;
   }
}

void initSpatialCellCoordinates(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   vector<CellID> cells = mpiGrid.get_cells();
   #pragma omp parallel for
   for (size_t i=0; i<cells.size(); ++i) {
      std::array<double, 3> cell_min = mpiGrid.geometry.get_min(cells[i]);
      std::array<double, 3> cell_length = mpiGrid.geometry.get_length(cells[i]);

      mpiGrid[cells[i]]->parameters[CellParams::XCRD] = cell_min[0];
      mpiGrid[cells[i]]->parameters[CellParams::YCRD] = cell_min[1];
      mpiGrid[cells[i]]->parameters[CellParams::ZCRD] = cell_min[2];
      mpiGrid[cells[i]]->parameters[CellParams::DX  ] = cell_length[0];
      mpiGrid[cells[i]]->parameters[CellParams::DY  ] = cell_length[1];
      mpiGrid[cells[i]]->parameters[CellParams::DZ  ] = cell_length[2];

      mpiGrid[cells[i]]->parameters[CellParams::CELLID] = cells[i];
      mpiGrid[cells[i]]->parameters[CellParams::REFINEMENT_LEVEL] = mpiGrid.get_refinement_level(cells[i]);
   }
}

/*
Record for each cell which processes own one or more of its face neighbors
 */
void setFaceNeighborRanks( dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid ) {

   const vector<CellID>& cells = getLocalCells();
   // TODO: Try a #pragma omp parallel for
   for (const auto& cellid : cells) {

      if (cellid == INVALID_CELLID) continue;

      SpatialCell* cell = mpiGrid[cellid];

      if (!cell) continue;

      cell->face_neighbor_ranks.clear();

      for (const auto& [neighbor, dir] : mpiGrid.get_face_neighbors_of(cellid)) {

         int neighborhood;

         // We store rank numbers into a map that has neighborhood ids as its key values.

         switch (dir) {
         case -3:
            neighborhood = SHIFT_M_Z_NEIGHBORHOOD_ID;
            break;
         case -2:
            neighborhood = SHIFT_M_Y_NEIGHBORHOOD_ID;
            break;
         case -1:
            neighborhood = SHIFT_M_X_NEIGHBORHOOD_ID;
            break;
         case +1:
            neighborhood = SHIFT_P_X_NEIGHBORHOOD_ID;
            break;
         case +2:
            neighborhood = SHIFT_P_Y_NEIGHBORHOOD_ID;
            break;
         case +3:
            neighborhood = SHIFT_P_Z_NEIGHBORHOOD_ID;
            break;
         default:
            cerr << "Invalid face neighbor dimension: " << dir << " in " << __FILE__ << ":" << __LINE__ << std::endl;
            abort();
         }

         cell->face_neighbor_ranks[neighborhood].insert(mpiGrid.get_process(neighbor));

      }
   }
}

void balanceLoad(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, SysBoundary& sysBoundaries){
   // Invalidate cached cell lists
   Parameters::meshRepartitioned = true;

   // tell other processes which velocity blocks exist in remote spatial cells
   phiprof::Timer balanceLoadTimer {"Balancing load", {"Load balance"}};

   // memory deallocation for dAMR. Do we want this got GPUs as well?
   phiprof::Timer deallocTimer {"deallocate boundary data"};
   //deallocate blocks in remote cells to decrease memory load
   deallocateRemoteCellBlocks(mpiGrid);
   deallocTimer.stop();

   //set weights based on each cells LB weight counter
   const vector<CellID>& cells = getLocalCells();
   for (size_t i=0; i<cells.size(); ++i){
      // Set cell weight. We could use different counters or number of blocks if different solvers are active.
      // if (P::propagateVlasovAcceleration)
      // When using the FS-SPLIT functionality, Jaro Hokkanen reported issues with using the regular
      // CellParams::LBWEIGHTCOUNTER, so use of blockscounts + 1 might be required.
      mpiGrid.set_cell_weight(cells[i], (Real)1 + mpiGrid[cells[i]]->parameters[CellParams::LBWEIGHTCOUNTER]);
   }

   phiprof::Timer initLBTimer {"dccrg.initialize_balance_load"};
   mpiGrid.initialize_balance_load(true);
   initLBTimer.stop();

   const std::unordered_set<CellID>& incoming_cells = mpiGrid.get_cells_added_by_balance_load();
   std::vector<CellID> incoming_cells_list (incoming_cells.begin(),incoming_cells.end());

   const std::unordered_set<CellID>& outgoing_cells = mpiGrid.get_cells_removed_by_balance_load();
   std::vector<CellID> outgoing_cells_list (outgoing_cells.begin(),outgoing_cells.end());

   /*transfer cells in parts to preserve memory*/
   phiprof::Timer transfersTimer {"Data transfers"};
   const uint64_t num_part_transfers=5;
   for (uint64_t transfer_part=0; transfer_part<num_part_transfers; transfer_part++) {
      //Set transfers on/off for the incoming cells in this transfer set and prepare for receive
      for (unsigned int i=0;i<incoming_cells_list.size();i++){
         CellID cell_id=incoming_cells_list[i];
         SpatialCell* cell = mpiGrid[cell_id];
         if (cell_id%num_part_transfers!=transfer_part) {
            cell->set_mpi_transfer_enabled(false);
         } else {
            cell->set_mpi_transfer_enabled(true);
         }
      }

      //Set transfers on/off for the outgoing cells in this transfer set
      for (unsigned int i=0; i<outgoing_cells_list.size(); i++) {
         CellID cell_id=outgoing_cells_list[i];
         SpatialCell* cell = mpiGrid[cell_id];
         if (cell_id%num_part_transfers!=transfer_part) {
            cell->set_mpi_transfer_enabled(false);
         } else {
            cell->set_mpi_transfer_enabled(true);
         }
      }

      for (size_t popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         // Set active population
         SpatialCell::setCommunicatedSpecies(popID);

         //Transfer velocity block list
         SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE1);
         mpiGrid.continue_balance_load();
         SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE2);
         mpiGrid.continue_balance_load();

         int prepareReceives {phiprof::initializeTimer("Preparing receives")};
         int receives = 0;
         #pragma omp parallel for schedule(guided)
         for (unsigned int i=0; i<incoming_cells_list.size(); i++) {
            CellID cell_id=incoming_cells_list[i];
            SpatialCell* cell = mpiGrid[cell_id];
            if (cell_id % num_part_transfers == transfer_part) {
               receives++;
               // reserve space for velocity block data in arriving remote cells
               phiprof::Timer timer {prepareReceives};
               cell->setNewSizeClear(popID);
               cell->prepare_to_receive_blocks(popID);
               timer.stop(1, "Spatial cells");
            }
         }
         if(receives == 0) {
            //empty phiprof timer, to avoid unneccessary divergence in unique
            //profiles (keep order same)
            phiprof::Timer timer {prepareReceives};
            timer.stop(0, "Spatial cells");
         }

         //do the actual transfer of data for the set of cells to be transferred
         phiprof::Timer transferTimer {"transfer_all_data"};
         SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
         mpiGrid.continue_balance_load();
         transferTimer.stop();

         // Free memory for cells that have been sent (the block data)
         for (unsigned int i=0;i<outgoing_cells_list.size();i++){
            CellID cell_id=outgoing_cells_list[i];
            SpatialCell* cell = mpiGrid[cell_id];

            // Free memory of this cell as it has already been transferred,
            // it will not be used anymore. NOTE: Only clears memory allocated
            // to the active population.
            if (cell_id % num_part_transfers == transfer_part) cell->clear(popID);
         }
      } // for-loop over populations
   } // for-loop over transfer parts
   transfersTimer.stop();

   //finish up load balancing
   phiprof::Timer finishLBTimer {"dccrg.finish_balance_load"};
   mpiGrid.finish_balance_load();
   finishLBTimer.stop();

   //Make sure transfers are enabled for all cells
   recalculateLocalCellsCache();
   #pragma omp parallel for
   for (uint i=0; i<cells.size(); ++i) {
      mpiGrid[cells[i]]->set_mpi_transfer_enabled(true);
   }

   // flag transfers if AMR
   phiprof::Timer computeTransferTimer {"compute_amr_transfer_flags"};
   flagSpatialCellsForAmrCommunication(mpiGrid,cells);
   computeTransferTimer.stop();

   // Communicate all spatial data for FULL neighborhood, which
   // includes all data with the exception of dist function data
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_copies_of_remote_neighbors(FULL_NEIGHBORHOOD_ID);

   phiprof::Timer updateBlocksTimer {"update block lists"};
   //new partition, re/initialize blocklists of remote cells.
   for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      updateRemoteVelocityBlockLists(mpiGrid,popID);
   }
   updateBlocksTimer.stop();

   phiprof::Timer updateBoundariesTimer {"update sysboundaries"};
   sysBoundaries.updateSysBoundariesAfterLoadBalance( mpiGrid );
   updateBoundariesTimer.stop();

   phiprof::Timer initSolversTimer {"Init solvers"};
   // Initialize field propagator (only if in use):
   if (Parameters::propagateField == true) {
      if (initializeFieldPropagatorAfterRebalance() == false) {
         logFile << "(MAIN): Field propagator did not initialize correctly!" << endl << writeVerbose;
         exit(1);
      }
   }
   initSolversTimer.stop();

   // Record ranks of face neighbors
   phiprof::Timer setRanksTimer("set face neighbor ranks");
   setFaceNeighborRanks( mpiGrid );
   setRanksTimer.stop();

   // Prepare cellIDs and pencils for AMR translation
   phiprof::Timer seedIdsTimer("GetSeedIdsAndBuildPencils");
   prepareSeedIdsAndPencils(mpiGrid);
   seedIdsTimer.stop();

#ifdef USE_GPU
   phiprof::Timer gpuMallocTimer("GPU_malloc");
   uint gpuMaxBlockCount = 0;
   vmesh::LocalID gpuBlockCount = 0;
   // Not parallelized
   const vector<CellID>& newCells = getLocalCells();
   const std::vector<CellID>& remote_cells = mpiGrid.get_remote_cells_on_process_boundary(FULL_NEIGHBORHOOD_ID);
   for (uint i=0; i<newCells.size()+remote_cells.size(); ++i) {
      SpatialCell* SC;
      if (i < newCells.size()) {
         SC = mpiGrid[newCells[i]];
      } else {
         SC = mpiGrid[remote_cells[i - newCells.size()]];
      }
      for (size_t popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
         const vmesh::VelocityMesh* vmesh = SC->get_velocity_mesh(popID);
         vmesh::VelocityBlockContainer* blockContainer = SC->get_velocity_blocks(popID);
         gpuBlockCount = vmesh->size();
         // checks if increased allocation is necessary, also performs deallocation first if necessary
         blockContainer->setNewCapacity(gpuBlockCount);
         if (gpuBlockCount > gpuMaxBlockCount) {
            gpuMaxBlockCount = gpuBlockCount;
         }
         // Ensure cell has sufficient reservation, then apply it
         SC->setReservation(popID,gpuBlockCount);
         SC->applyReservation(popID);
         SC->dev_upload_population(popID);
      }
   }
   // Call GPU routines for per-thread memory allocation for Vlasov solvers
   // deallocates first if necessary
   //GPUTODO: Also count how many pencils exist
   gpu_vlasov_allocate(gpuMaxBlockCount);
   gpu_acc_allocate(gpuMaxBlockCount);
   gpuMallocTimer.stop();
#endif // end USE_GPU
}

/*
  Adjust sparse velocity space to make it consistent in all 6 dimensions.

  Further documentation in grid.h
*/
bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const vector<CellID>& cellsToAdjust,
                          bool doPrepareToReceiveBlocks,
                          const uint popID) {
   phiprof::Timer readjustBlocksTimer {"re-adjust blocks", {"Block adjustment"}};
   SpatialCell::setCommunicatedSpecies(popID);

   const vector<CellID>& cells = getLocalCells();
   // Batch call
   update_velocity_block_content_lists(mpiGrid,cells,popID);

   if (doPrepareToReceiveBlocks) {
      // We are in the last substep of acceleration, so need to account for neighbours
      phiprof::Timer transferTimer {"Transfer with_content_list", {"MPI"}};
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_WITH_CONTENT_STAGE1 );
      mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_WITH_CONTENT_STAGE2 );
      mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);
      transferTimer.stop();
   }
   // Batch adjusts velocity blocks in local spatial cells, doesn't adjust velocity blocks in remote cells.
   adjust_velocity_blocks_in_cells(mpiGrid, cellsToAdjust, popID, doPrepareToReceiveBlocks);

   // prepare to receive block data
   if (doPrepareToReceiveBlocks) {
      updateRemoteVelocityBlockLists(mpiGrid,popID);
   }
   return true;
}

/*! Shrink to fit velocity space data to save memory.
 * \param mpiGrid Spatial grid
 */
void shrink_to_fit_grid_data(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   const std::vector<CellID>& cells = getLocalCells();
   const std::vector<CellID>& remote_cells = mpiGrid.get_remote_cells_on_process_boundary(FULL_NEIGHBORHOOD_ID);
   #pragma omp parallel for
   for(size_t i=0; i<cells.size() + remote_cells.size(); ++i) {
      if(i < cells.size()){
         SpatialCell* target= mpiGrid[cells[i]];
         if (target!=nullptr){
            target->shrink_to_fit();
         }
      }else{
         SpatialCell* target= mpiGrid[remote_cells[i - cells.size()]];
         if (target!=nullptr){
            target->shrink_to_fit();
         }
      }
   }
}

/*! Estimates memory consumption and writes it into logfile. Collective operation on MPI_COMM_WORLD
 * \param mpiGrid Spatial grid
 */
void report_grid_memory_consumption(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   /*now report memory consumption into logfile*/
   const vector<CellID>& cells = getLocalCells();
   const std::vector<CellID> remote_cells = mpiGrid.get_remote_cells_on_process_boundary();
   int rank,n_procs;
   MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   /* Compute memory statistics of the memory consumption of the spatial cells.
    * Internally we use double as MPI does
    * not define proper uint64_t datatypes for MAXLOCNot Real, as we
    * want double here not to loose accuracy.
    */

   /*report data for memory needed by blocks*/
   double mem[6] = {0};
   double sum_mem[6];

   for(unsigned int i=0;i<cells.size();i++){
      mem[0] += mpiGrid[cells[i]]->get_cell_memory_size();
      mem[3] += mpiGrid[cells[i]]->get_cell_memory_capacity();
   }

   for(unsigned int i=0;i<remote_cells.size();i++){
      mem[1] += mpiGrid[remote_cells[i]]->get_cell_memory_size();
      mem[4] += mpiGrid[remote_cells[i]]->get_cell_memory_capacity();
   }

   mem[2] = mem[0] + mem[1];//total meory according to size()
   mem[5] = mem[3] + mem[4];//total memory according to capacity()


   MPI_Reduce(mem, sum_mem, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   logFile << "(MEM) Total size: " << sum_mem[2] << endl;
   logFile << "(MEM) Total capacity " << sum_mem[5] << endl;

   struct {
      double val;
      int   rank;
   } max_mem[3],mem_usage_loc[3],min_mem[3];
   for(uint i = 0; i<3; i++){
      mem_usage_loc[i].val = mem[i + 3]; //report on capacity numbers (6: local cells, 7: remote cells, 8: all cells)
      mem_usage_loc[i].rank = rank;
   }

   MPI_Reduce(mem_usage_loc, max_mem, 3, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   MPI_Reduce(mem_usage_loc, min_mem, 3, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

   logFile << "(MEM)   Average capacity: " << sum_mem[5]/n_procs << " local cells " << sum_mem[3]/n_procs << " remote cells " << sum_mem[4]/n_procs << endl;
   logFile << "(MEM)   Max capacity:     " << max_mem[2].val   << " on  process " << max_mem[2].rank << endl;
   logFile << "(MEM)   Min capacity:     " << min_mem[2].val   << " on  process " << min_mem[2].rank << endl;
   logFile << writeVerbose;
}

/*! Deallocates all block data in remote cells in order to save
 *  memory
 * \param mpiGrid Spatial grid
 */
void deallocateRemoteCellBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid) {
   const std::vector<uint64_t> incoming_cells
      = mpiGrid.get_remote_cells_on_process_boundary(VLASOV_SOLVER_NEIGHBORHOOD_ID);
   for(unsigned int i=0;i<incoming_cells.size();i++){
      uint64_t cell_id=incoming_cells[i];
      SpatialCell* cell = mpiGrid[cell_id];
      if (cell != NULL) {
         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID)
            cell->clear(popID);
      }
   }

}

/*
Updates velocity block lists between remote neighbors and prepares local
copies of remote neighbors for receiving velocity block data.
*/
void updateRemoteVelocityBlockLists(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const uint popID,
   const uint neighborhood/*=DIST_FUNC_NEIGHBORHOOD_ID default*/
)
{
   SpatialCell::setCommunicatedSpecies(popID);

   // update velocity block lists For small velocity spaces it is
   // faster to do it in one operation, and not by first sending size,
   // then list. For large we do it in two steps
   phiprof::Timer updateTimer {"Velocity block list update", {"MPI"}};
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE1);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);
   SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE2);
   mpiGrid.update_copies_of_remote_neighbors(neighborhood);
   updateTimer.stop();

   // Prepare spatial cells for receiving velocity block data
   phiprof::Timer receivesTimer {"Preparing receives"};
   const std::vector<uint64_t> incoming_cells = mpiGrid.get_remote_cells_on_process_boundary(neighborhood);

   #pragma omp parallel for
   for (unsigned int i=0; i<incoming_cells.size(); ++i) {
     uint64_t cell_id = incoming_cells[i];
     SpatialCell* cell = mpiGrid[cell_id];
     if (cell == NULL) {
        #ifdef DEBUG_VLASIATOR
        for (const auto& cell: mpiGrid.local_cells) {
           if (cell.id == cell_id) {
              cerr << __FILE__ << ":" << __LINE__ << std::endl;
              abort();
           }
           for (const auto& neighbor: cell.neighbors_of) {
              if (neighbor.id == cell_id) {
                 cerr << __FILE__ << ":" << __LINE__ << std::endl;
                 abort();
              }
           }
        }
        #endif
        continue;
     }
     cell->setNewSizeClear(popID);
     cell->prepare_to_receive_blocks(popID);
   }

   receivesTimer.stop(incoming_cells.size(), "SpatialCells");
}

/*
  Set stencils. These are the stencils (in 2D, real ones in 3D of
  course). x are stencil neighbor to cell local cell o:

NEAREST SYSBOUNDARIES  (nearest neighbor)
-----------
  xxx
  xox
  xxx
-----------

EXTENDED_SYSBOUNDARIES (second nearest neighbor, also in diagonal)
-----------
  xxxxx
  xxxxx
  xxoxx
  xxxxx
  xxxxx
-----------

VLASOV
-----------
    x
    x
  xxoxx
    x
    x
-----------

VLASOV_{XYZ}
-----------
 xxoxxx
-----------

VLASOV_TARGET_{XYZ}
-----------
  xox

-----------

DIST_FUNC  (Includes all cells which should know about each others blocks and have space for them. VLASOV + SYSBOUNDARIES.
-----------
    x
   xxx
  xxoxx
   xxx
    x

-----------


FULL (Includes all possible communication, possible AMR extension)
-----------
    A
  xxxxx
  xxxxx
 AxxoxxA
  xxxxx
  xxxxx
    A
-----------

SHIFT_M_X    ox
SHIFT_P_X   xo
 Y, Z in the same way
*/

void initializeStencils(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid){
   // set reduced neighborhoods
   typedef dccrg::Types<3>::neighborhood_item_t neigh_t;

   // set a reduced neighborhood for nearest neighbours
   std::vector<neigh_t> neighborhood;
   for (int z = -1; z <= 1; z++) {
      for (int y = -1; y <= 1; y++) {
         for (int x = -1; x <= 1; x++) {
            if (x == 0 && y == 0 && z == 0) {
               continue;
            }
            neigh_t offsets = {{x, y, z}};
            neighborhood.push_back(offsets);
         }
      }
   }
   //mpiGrid.add_neighborhood(FIELD_SOLVER_NEIGHBORHOOD_ID, neighborhood);
   if (!mpiGrid.add_neighborhood(NEAREST_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood NEAREST_NEIGHBORHOOD_ID \n";
      abort();
   }
   if (!mpiGrid.add_neighborhood(SYSBOUNDARIES_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SYSBOUNDARIES_NEIGHBORHOOD_ID \n";
      abort();
   }

   neighborhood.clear();
   for (int z = -2; z <= 2; z++) {
      for (int y = -2; y <= 2; y++) {
         for (int x = -2; x <= 2; x++) {
            if (x == 0 && y == 0 && z == 0) {
               continue;
            }
            neigh_t offsets = {{x, y, z}};
            neighborhood.push_back(offsets);
         }
      }
   }
   if (!mpiGrid.add_neighborhood(SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SYSBOUNDARIES_EXTENDED_NEIGHBORHOOD_ID \n";
      abort();
   }

   // In spatial AMR using DCCRG, the neighbors are considered relative to a given cell's size.
   // To get two coarse neighbors from a fine cell at interfaces, the stencil size needs to be increased by one.
   int addStencilDepth = 0;
   if (P::amrMaxSpatialRefLevel > 0) {
      switch (VLASOV_STENCIL_WIDTH) {
	 case 1:
	    // Required cells will be included already
	    break;
	 case 2:
	    // looking from high to low refinement: stencil 2 will only give 1 cell, so need to add 1
	    addStencilDepth = 1;
            break;
         case 3:
	    // looking from high to low refinement: stencil 3 will only give 2 cells, so need to add 2
	    // to reach surely into the third low-refinement neighbour
            addStencilDepth = 2;
            break;
         default:
            std::cerr<<"Warning: unrecognized VLASOV_STENCIL_WIDTH in grid.cpp"<<std::endl;
      }
   }
   int full_neighborhood_size = max(2, VLASOV_STENCIL_WIDTH);

   neighborhood.clear();
   for (int z = -full_neighborhood_size; z <= full_neighborhood_size; z++) {
      for (int y = -full_neighborhood_size; y <= full_neighborhood_size; y++) {
         for (int x = -full_neighborhood_size; x <= full_neighborhood_size; x++) {
            if (x == 0 && y == 0 && z == 0) {
               continue;
            }
            neigh_t offsets = {{x, y, z}};
            neighborhood.push_back(offsets);
         }
      }
   }
   /* Add extra face neighbors if required by AMR */
   for (int d = full_neighborhood_size+1; d <= full_neighborhood_size+addStencilDepth; d++) {
      neighborhood.push_back({{ d, 0, 0}});
      neighborhood.push_back({{-d, 0, 0}});
      neighborhood.push_back({{0, d, 0}});
      neighborhood.push_back({{0,-d, 0}});
      neighborhood.push_back({{0, 0, d}});
      neighborhood.push_back({{0, 0,-d}});
   }
   /*all possible communication pairs*/
   if( !mpiGrid.add_neighborhood(FULL_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood FULL_NEIGHBORHOOD_ID \n";
      abort();
   }

   /*stencils for semilagrangian propagators*/
   neighborhood.clear();
   for (int d = -VLASOV_STENCIL_WIDTH-addStencilDepth; d <= VLASOV_STENCIL_WIDTH+addStencilDepth; d++) {
     if (d != 0) {
        neighborhood.push_back({{d, 0, 0}});
        neighborhood.push_back({{0, d, 0}});
        neighborhood.push_back({{0, 0, d}});
     }
   }
   if( !mpiGrid.add_neighborhood(VLASOV_SOLVER_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood VLASOV_SOLVER_NEIGHBORHOOD_ID \n";
      abort();
   }

   // add remaining nearest neighbors for DIST_FUNC neighborhood
   for (int z = -1; z <= 1; z++) {
      for (int y = -1; y <= 1; y++) {
         for (int x = -1; x <= 1; x++) {
            //do not add cells already in neighborhood (vlasov solver)
            if (x == 0 && y == 0 ) continue;
            if (x == 0 && z == 0 ) continue;
            if (y == 0 && z == 0 ) continue;

            neigh_t offsets = {{x, y, z}};
            neighborhood.push_back(offsets);
         }
      }
   }
   if( !mpiGrid.add_neighborhood(DIST_FUNC_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood DIST_FUNC_NEIGHBORHOOD_ID \n";
      abort();
   }

   neighborhood.clear();
   for (int d = -VLASOV_STENCIL_WIDTH-addStencilDepth; d <= VLASOV_STENCIL_WIDTH+addStencilDepth; d++) {
     if (d != 0) {
        neighborhood.push_back({{d, 0, 0}});
     }
   }
   if( !mpiGrid.add_neighborhood(VLASOV_SOLVER_X_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood VLASOV_SOLVER_X_NEIGHBORHOOD_ID \n";
      abort();
   }


   neighborhood.clear();
   for (int d = -VLASOV_STENCIL_WIDTH-addStencilDepth; d <= VLASOV_STENCIL_WIDTH+addStencilDepth; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, d, 0}});
     }
   }
   if( !mpiGrid.add_neighborhood(VLASOV_SOLVER_Y_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood VLASOV_SOLVER_Y_NEIGHBORHOOD_ID \n";
      abort();
   }


   neighborhood.clear();
   for (int d = -VLASOV_STENCIL_WIDTH-addStencilDepth; d <= VLASOV_STENCIL_WIDTH+addStencilDepth; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, 0, d}});
     }
   }
   if (!mpiGrid.add_neighborhood(VLASOV_SOLVER_Z_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood VLASOV_SOLVER_Z_NEIGHBORHOOD_ID \n";
      abort();
   }

   neighborhood.clear();
   for (int d = -1; d <= 1; d++) {
     if (d != 0) {
        neighborhood.push_back({{d, 0, 0}});
     }
   }
   if (!mpiGrid.add_neighborhood(VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood VLASOV_SOLVER_TARGET_X_NEIGHBORHOOD_ID \n";
      abort();
   }

   neighborhood.clear();
   for (int d = -1; d <= 1; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, d, 0}});
     }
   }
   if (!mpiGrid.add_neighborhood(VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood VLASOV_SOLVER_TARGET_Y_NEIGHBORHOOD_ID \n";
      abort();
   }

   neighborhood.clear();
   for (int d = -1; d <= 1; d++) {
     if (d != 0) {
        neighborhood.push_back({{0, 0, d}});
     }
   }
   if (!mpiGrid.add_neighborhood(VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood VLASOV_SOLVER_TARGET_Z_NEIGHBORHOOD_ID \n";
      abort();
   }


   neighborhood.clear();
   neighborhood.push_back({{1, 0, 0}});
   if (!mpiGrid.add_neighborhood(SHIFT_M_X_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SHIFT_M_X_NEIGHBORHOOD_ID \n";
      abort();
   }
   neighborhood.clear();
   neighborhood.push_back({{0, 1, 0}});
   if (!mpiGrid.add_neighborhood(SHIFT_M_Y_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SHIFT_M_Y_NEIGHBORHOOD_ID \n";
      abort();
   }
   neighborhood.clear();
   neighborhood.push_back({{0, 0, 1}});
   if (!mpiGrid.add_neighborhood(SHIFT_M_Z_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SHIFT_M_Z_NEIGHBORHOOD_ID \n";
      abort();
   }
   neighborhood.clear();
   neighborhood.push_back({{-1, 0, 0}});
   if (!mpiGrid.add_neighborhood(SHIFT_P_X_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SHIFT_P_X_NEIGHBORHOOD_ID \n";
      abort();
   }
   neighborhood.clear();
   neighborhood.push_back({{0, -1, 0}});
   if (!mpiGrid.add_neighborhood(SHIFT_P_Y_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SHIFT_P_Y_NEIGHBORHOOD_ID \n";
      abort();
   }
   neighborhood.clear();
   neighborhood.push_back({{0, 0, -1}});
   if (!mpiGrid.add_neighborhood(SHIFT_P_Z_NEIGHBORHOOD_ID, neighborhood)){
      std::cerr << "Failed to add neighborhood SHIFT_P_Z_NEIGHBORHOOD_ID \n";
      abort();
   }
}

void mapRefinement(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid) {
   phiprof::Timer timer {"Map Refinement Level to FsGrid"};
   const FsGridTools::FsIndex_t *localDims = &technicalGrid.getLocalSize()[0];

   // #pragma omp parallel for collapse(3)
   for (FsGridTools::FsIndex_t k=0; k<localDims[2]; k++) {
      for (FsGridTools::FsIndex_t j=0; j<localDims[1]; j++) {
         for (FsGridTools::FsIndex_t i=0; i<localDims[0]; i++) {

            const std::array<FsGridTools::FsIndex_t, 3> mapIndices = technicalGrid.getGlobalIndices(i,j,k);
            const dccrg::Types<3>::indices_t  indices = {{(uint64_t)mapIndices[0],(uint64_t)mapIndices[1],(uint64_t)mapIndices[2]}}; //cast to avoid warnings
            CellID dccrgCellID2 = mpiGrid.get_existing_cell(indices, 0, mpiGrid.mapping.get_maximum_refinement_level());
            int amrLevel= mpiGrid.get_refinement_level(dccrgCellID2);
            technicalGrid.get(i, j, k)-> refLevel =amrLevel ;
         }
      }
   }
}

bool adaptRefinement(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> & technicalGrid, SysBoundary& sysBoundaries, Project& project, int useStatic) {
   phiprof::Timer amrTimer {"Re-refine spatial cells"};
   int refines {0};
   if (useStatic > -1) {
      project.forceRefinement(mpiGrid, useStatic);
   } else {
      // Restarts don't have all the data needed to calculate indices so they are read directly
      if (P::tstep != P::tstep_min) {
         calculateScaledDeltasSimple(mpiGrid);
      }
      SpatialCell::set_mpi_transfer_type(Transfer::REFINEMENT_PARAMETERS);
      mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);

      refines = project.adaptRefinement(mpiGrid);
   }

   int cells = getLocalCells().size();
   MPI_Allreduce(MPI_IN_PLACE, &refines, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   double ratio_refines = static_cast<double>(refines) / static_cast<double>(cells);
   logFile << "(AMR) Refining " << refines << " cells, " << 100.0 * ratio_refines << "% of grid" << std::endl;

   phiprof::Timer dccrgTimer {"dccrg refinement"};

   phiprof::Timer initTimer {"initialize refines"};
   mpiGrid.initialize_refines();
   initTimer.stop();

   refines = mpiGrid.get_local_cells_to_refine().size();
   double coarsens = mpiGrid.get_local_cells_to_unrefine().size();
   MPI_Allreduce(MPI_IN_PLACE, &refines, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &coarsens, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   ratio_refines = static_cast<double>(refines) / static_cast<double>(cells);
   double ratio_coarsens = static_cast<double>(coarsens) / static_cast<double>(cells);
   logFile << "(AMR) Refining " << refines << " cells after induces, " << 100.0 * ratio_refines << "% of grid" << std::endl;
   logFile << "(AMR) Coarsening " << coarsens << " cells after induces, " << 100.0 * ratio_coarsens << "% of grid" << std::endl;

   double newBytes{0};
   phiprof::Timer estimateMemoryTimer {"Estimate memory usage"};
   for (auto id : mpiGrid.get_local_cells_to_refine()) {
      newBytes += 8 * mpiGrid[id]->get_cell_memory_capacity();
   }

   // Rougher estimate than above
   // Unrefined cells have a transitive memory footprint since parent and children exist at same time
   for (auto id : mpiGrid.get_local_cells_to_unrefine()) {
      newBytes += mpiGrid[id]->get_cell_memory_capacity();
   }

   report_process_memory_consumption(newBytes);
   estimateMemoryTimer.stop();

   logFile.flush(false);

   // Bailout from estimate
   // clunky...
   int bailout {0};
   phiprof::Timer bailoutAllreduceTimer {"Bailout-allreduce"};
   MPI_Allreduce(&(globalflags::bailingOut), &bailout, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   bailoutAllreduceTimer.stop();

   if (bailout) {
      return false;
   }

   // New cells created by refinement
   phiprof::Timer executeTimer {"execute refines"};
   auto newChildren = mpiGrid.execute_refines();
   executeTimer.stop();

   std::vector<CellID> receives;
   for (auto const& [key, val] : mpiGrid.get_cells_to_receive()) {
      for (auto i : val) {
         receives.push_back(i.first);
      }
   }

   phiprof::Timer transfersTimer {"transfers"};
   for (size_t popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      // Set active population
      SpatialCell::setCommunicatedSpecies(popID);

      //Transfer velocity block list
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE1);
      mpiGrid.continue_refining();
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_LIST_STAGE2);
      mpiGrid.continue_refining();

      int prepareReceives {phiprof::initializeTimer("Preparing receives")};
      for (CellID id : receives) {
         // reserve space for velocity block data in arriving remote cells
         phiprof::Timer timer {prepareReceives};
         mpiGrid[id]->setNewSizeClear(popID);
         mpiGrid[id]->prepare_to_receive_blocks(popID);
         timer.stop(1, "Spatial cells");
      }

      if(receives.empty()) {
         //empty phiprof timer, to avoid unneccessary divergence in unique
         //profiles (keep order same)
         phiprof::Timer timer {prepareReceives};
         timer.stop(0, "Spatial cells");
      }

      //do the actual transfer of data for the set of cells to be transferred
      phiprof::Timer transferTimer {"transfer_all_data"};
      SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
      mpiGrid.continue_refining();
      transferTimer.stop();
   }
   transfersTimer.stop();

   phiprof::Timer copyChildrenTimer {"copy to children"};
   for (CellID id : newChildren) {
      *mpiGrid[id] = *mpiGrid[mpiGrid.get_parent(id)];
      // Irrelevant?
      mpiGrid[id]->parameters[CellParams::AMR_ALPHA1] /= 2.0;
      mpiGrid[id]->parameters[CellParams::AMR_ALPHA2] /= 2.0;
      mpiGrid[id]->parameters[CellParams::RECENTLY_REFINED] = 1;
   }
   copyChildrenTimer.stop(newChildren.size(), "Spatial cells");

   // Old cells removed by refinement
   phiprof::Timer copyParentsTimer {"copy to parents"};
   std::set<CellID> processed;
   for (CellID id : mpiGrid.get_removed_cells()) {
      CellID parent = mpiGrid.get_existing_cell(mpiGrid.get_center(id));
      if (!processed.count(parent)) {
         std::vector<CellID> children = mpiGrid.get_all_children(parent);
         // Make sure cell contents aren't garbage
         *mpiGrid[parent] = *mpiGrid[id];

         for (uint popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
            SBC::averageCellData(mpiGrid, children, mpiGrid[parent], popID, 1);
         }

         // Averaging moments
         calculateCellMoments(mpiGrid[parent], true, false);

         processed.insert(parent);
      }
   }
   copyParentsTimer.stop(processed.size(), "Spatial cells");

   phiprof::Timer finishTimer {"finish refining"};
   mpiGrid.finish_refining();
   finishTimer.stop();
   dccrgTimer.stop();

   recalculateLocalCellsCache();
   initSpatialCellCoordinates(mpiGrid);

   SpatialCell::set_mpi_transfer_type(Transfer::CELL_DIMENSIONS);
   mpiGrid.update_copies_of_remote_neighbors(SYSBOUNDARIES_NEIGHBORHOOD_ID);

   mapRefinement(mpiGrid, technicalGrid);

   // Initialise system boundary conditions (they need the initialised positions!!)
	// This needs to be done before LB
   sysBoundaries.classifyCells(mpiGrid,technicalGrid);

   //SpatialCell::set_mpi_transfer_type(Transfer::ALL_DATA);
   SpatialCell::set_mpi_transfer_type(Transfer::CELL_PARAMETERS);
   mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);

   // Is this needed?
   technicalGrid.updateGhostCells(); // This needs to be done at some point
   for (size_t popID=0; popID<getObjectWrapper().particleSpecies.size(); ++popID) {
      updateRemoteVelocityBlockLists(mpiGrid, popID, NEAREST_NEIGHBORHOOD_ID);
   }

   if (P::shouldFilter) {
      project.filterRefined(mpiGrid);
   }
   return true;
}
