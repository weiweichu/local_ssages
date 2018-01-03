/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Weiwei Chu <weiweichu@uchicago.edu>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "Meta.h"

namespace SSAGES
{
  //! Well-tempered Metadynamics sampling method
  /*!
   * \ingroup Methods
   * The notation used here is drawn largely fromA. Barducci, G. Bussi and M.
   * Parrinello, Phys. Rev. Lett. 100, 020603 (2008).
   */
  class WellTemperedMeta : public Meta
  {
    private:
      double delta_T_ = 0.0;

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Current simulation snapshot.
		 * \param cvmanager Collective variable manager.
		 */
		  void PreSimulation(Snapshot* snapshot, const class CVManager& cvmanager) override;

      //! Adds a new hill.
      /*!
       * \param cvs List of CVs.
       * \param iteration Current iteration.
       */
      void AddHill(const CVList& cvs, int iteration) override;

      //! Computes the bias force.
       /*!
        * \param cvs List of CVs.
        */
      void CalcBiasForce(const CVList& cvs) override;

    public:
      //! Constructor
      /*!
       * \param world MPI global communicator.
       * \param comm MPI local communicator.
       * \param height Initial height of the hills to be deposited.
       * \param widths Width of the hills to be deposited.
       * \param hillfreq Frequency of depositing hills.
       * \param frequency Frequency of invoking this method.
       * \param delta_T Temperature offset when adjusting hill heights.
       *
       * Create instance of WellTemperedMeta method.
       */
      WellTemperedMeta(const MPI_Comm& world,
                       const MPI_Comm& comm,
                       double height,
                       const std::vector<double>& widths,
                       const std::vector<double>& lowerb,
                       const std::vector<double>& upperb,
                       const std::vector<double>& lowerk,
                       const std::vector<double>& upperk,
                       Grid<Vector>* grid,
                       Grid<double>* grid_potential,
                       unsigned int hillfreq,
                       unsigned int potential_freq,
                       unsigned int frequency,
                       double delta_T) :
      Meta(world, comm, height, widths, lowerb, upperb, lowerk, upperk, grid, grid_potential, hillfreq, potential_freq, frequency),
      delta_T_(delta_T)
      {} 

		//! Destructor.
		~WellTemperedMeta() {}
	};
}
			
