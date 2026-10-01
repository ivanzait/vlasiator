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
 * 2D GEM-like magnetic reconnection setup (Birn et al. 2001).
 *
 * Geometry: x is the reconnection (outflow) direction,
 *           z is the current-sheet normal direction,
 *           y is the out-of-plane (guide-field) direction (y_length = 1, periodic).
 *
 * Magnetic field:
 *   Bx = B0 * tanh(z/lambda)  +  flux-function perturbation
 *   By = Bg  (guide field, 0 for anti-parallel)
 *   Bz =                         flux-function perturbation
 *   psi(x,z) = Psi0 * cos(2*pi*x/Lx) * cos(pi*z/Lz)
 *   -> Bx_pert = -dpsi/dz,  Bz_pert = dpsi/dx
 *
 * Density / temperature:
 *   n(z) = n0 * sech^2(z/lambda) + n_bg    (n_bg ~ 0.2 * n0)
 *   Isotropic Maxwellian with temperature T.
 *   Pressure balance requires: n0 * kB * T = B0^2/(2*mu0).
 */

#ifndef RECONNECTION_H
#define RECONNECTION_H

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {

struct ReconnectionSpeciesParameters {
   Real T;           // Temperature [K]
   Real n;           // Peak number density [m^-3]
   Real n_bg_frac;   // Background density fraction (n_bg = n_bg_frac * n)
   Real nSpaceSamples;
   Real nVelocitySamples;
};

class Reconnection : public TriAxisSearch {
public:
   Reconnection();
   virtual ~Reconnection();

   virtual bool initialize(void) override;
   virtual void addParameters(void) override;
   virtual void getParameters(void) override;
   virtual void calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) override;
   virtual void setProjectBField(fsgrids::perbspan perb,
                                 fsgrids::bgbspan bgb,
                                 fsgrids::technicalspan technical, FieldSolverGrid& fsgrid) override;

   virtual Realf fillPhaseSpace(spatial_cell::SpatialCell* cell,
                                const uint popID,
                                const uint nRequested) const override;
   virtual Realf probePhaseSpace(spatial_cell::SpatialCell* cell,
                                 const uint popID,
                                 Real vx_in, Real vy_in, Real vz_in) const override;
   virtual std::vector<std::array<Real, 3>> getV0(
      creal x, creal y, creal z, const uint popID) const override;

   // Public: RP::add<T> binds these member variables by reference at parse time.
   Real lambda;        // Current-sheet half-thickness [m]
   Real Bx0;           // Lobe field strength [T]
   Real By0;           // Guide field [T]
   Real Bz0;           // Normal field (usually 0 initially) [T]
   Real Psi0;          // Perturbation amplitude [T*m]

   std::vector<ReconnectionSpeciesParameters> speciesParams;
   std::vector<ReconnectionSpeciesParameters*> speciesParamsRead;

protected:
   Real getDensity(creal z, const uint popID) const;

}; // class Reconnection
} // namespace projects

#endif // RECONNECTION_H
