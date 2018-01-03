/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Jonathan K. Whitmer <jwhitme1@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
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
#include "StandardMeta.h"
#include "CVs/CVManager.h"
#include "Snapshot.h"
#include "spline.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
using namespace Json;

namespace SSAGES
{
  // Pre-simulation hook.
	void StandardMeta::PreSimulation(Snapshot* snapshot, const CVManager& cvmanager)
	{
		auto cvs = cvmanager.GetCVs(cvmask_);
    if (!potential_initial_) {
  		if(world_.rank() == 0)
  		{
  			hillsout_.open("hills.out");
  			hillsout_ << "#Iteration "; 
  
  			for(size_t i = 0; i < cvs.size(); ++i)
  				hillsout_ << "center." << i << " ";
  	
  			for(size_t i = 0; i < cvs.size(); ++i)
  				hillsout_ << "sigma." << i << " ";
  		  
  			hillsout_ << "height.";
        
        hillsout_ << std::endl;
  				
  			hillsout_.close();
       
        if (grid_ != nullptr) {
          potential_out_.open("potential.out");
          potential_out_.close();
        }
  		}
    }
	  
  // Initialize grid to zero and grid_potential to initial heights.
	 	if(grid_ != nullptr)
		{
			Vector vec = Vector::Zero(cvs.size());
			std::fill(grid_->begin(), grid_->end(), vec);
		} 

    // There are three ways to initialize grid_potential:
    // 1. Read in from potential files at a certain timestep.
    // 2. Load hills and calculate potential from hills
    // 3. Initialize to zeros.

    // Initialize grid_potential to read in initial potential.
    if (potential_initial_) {
    }

    // Initialize grid_potential by adding up hills.
    else if(hills_.size() > 0) {
  		int n = cvs.size();
      std::vector<double> dx(n, 0.0), dxp(n, 0.0), df(n, 1.0);
      double dp = 1.0;
      std::ofstream hill_state;
      hill_state.open("hill_state.out");
      int nhill = 0;
      for (auto& hill: hills_) {
        hill_state << nhill << std::endl;
        nhill ++;
        auto ip = grid_potential_->begin();
        auto it = grid_->begin();
        for(; it != grid_->end(); ++it, ++ip)
        {
            auto& val = *it;
            auto& valp = *ip;
            auto coord = it.coordinates(), coordp = ip.coordinates();
            // Compute difference between grid point and current val.
            for(size_t i = 0; i < n; ++i)
            {
                dx[i] = coord[i] - hill.center[i];
                dxp[i] = coordp[i] - hill.center[i];
                if (dx[i] - dxp[i] > 0.00001) {
                  std::cerr << "dx and dxp not the same" << std::endl;
                }
                df[i] = 1.;
            }
            dp = 1.0;
            // Compute derivative.
            for(size_t i = 0; i < n; ++i)
            {
                dp *= gaussian(dxp[i], hill.width[i]);
                for(size_t j = 0; j < n; ++j)
                {
                // Need to reduce the dimension of grid_potential by size of cvs.
                    if(j != i)
                      df[i] *= gaussian(dx[j], hill.width[j]);
                    else
                      df[i] *= gaussianDerv(dx[j], hill.width[j]);
                }
            }
            // Add to grid.
            for(size_t i = 0; i < n; ++i){
                val[i] += height_ *  df[i];
            }
            // Add to grid_potential.
            valp += height_ * dp;
        }
      }
    }

    // Initialize grid_potential to zero.
    else{
      std::fill(grid_potential_->begin(), grid_potential_->end(), 0);
    }

		auto n = snapshot->GetTargetIterations();
		n = n ? n : 1e5; // Pre-allocate at least something.
	
		hills_.reserve(n+1);
		//widths_.reserve(n+1);
    //heights_.reserve(n+1);
		derivatives_.resize(cvs.size());
		tder_.resize(cvs.size());
		dx_.resize(cvs.size());
	}


	// Drop a new hill.
	void StandardMeta::AddHill(const CVList& cvs, int iteration)
	{
		int n = cvs.size();

		// Assume we have the same number of procs per walker.
		int nwalkers = world_.size()/comm_.size();

		// We need to exchange CV values across the walkers 
		// and to each proc on a walker.	
		std::vector<double> cvals(n*nwalkers, 0);

		if(comm_.rank() == 0)
		{
			for(auto i = 0, j = world_.rank()/comm_.size()*n; i < n; ++i,++j)
				cvals[j] = cvs[i]->GetValue();
		}

		// Reduce across all processors and add hills.
		MPI_Allreduce(MPI_IN_PLACE, cvals.data(), n*nwalkers, MPI_DOUBLE, MPI_SUM, world_);
		
		for(int i = 0; i < n*nwalkers; i += n)
		{
			std::vector<double> cval(cvals.begin() + i, cvals.begin() + i + n);
			hills_.emplace_back(cval, widths_, height_);
			
			// Write hill to file.
			if(world_.rank() == 0)
				PrintHill(hills_.back(), iteration);
		}

		// If grid is defined, add bias onto grid. 
		if(grid_ != nullptr)
		{
			std::vector<double> dx(n, 0.0), dxp(n, 1.0), df(n, 1.0);
      double dp = 1.0;
			auto& hill = hills_.back();
      std::vector<double> cval(n, 0.0);
      bool inbounds = true;
      for(size_t i = 0; i < n; ++i) {
        cval[i] = cvs[i]->GetValue();
        if(cval[i] < grid_->GetLower(i) || cval[i] > grid_->GetUpper(i))
          inbounds = false;
      }
      if(!inbounds) {
        if (comm_.rank() == 0) {
          std::cerr << "WellTemperedMeta: out of bounds ( ";
          for (auto& v: cval)
            std::cerr << v << " ";
          std::cerr << ")" << std::endl;
        }
      }
      auto ip = grid_potential_->begin();
			for(auto it = grid_->begin(); it != grid_->end(); ++it, ++ip)
			{
				auto& val = *it;
        auto& valp = *ip;
				auto coord = it.coordinates(), coordp = ip.coordinates();

				// Compute difference between grid point and current val. 
				for(size_t i = 0; i < n; ++i)
				{
					dx[i] = -cvs[i]->GetDifference(coord[i]);
          dxp[i] = -cvs[i]->GetDifference(coordp[i]);
					df[i] = 1.;
				}
        dp = 1.0;

				// Compute derivative.
				for(size_t i = 0; i < n; ++i)
				{
          dp *= gaussian(dxp[i], hill.width[i]);
					for(size_t j = 0; j < n; ++j)
					{
						if(j != i) 
							df[i] *= gaussian(dx[j], hill.width[j]);
						else
							df[i] *= gaussianDerv(dx[j], hill.width[j]);
					}
				}

				// Add to grid. 
				for(size_t i = 0; i < n; ++i){
					val[i] += height_*df[i];
        }
        valp += height_ * dp;
			}
		}
	}


	void StandardMeta::CalcBiasForce(const CVList& cvs)
	{	
		// Reset bias and derivatives.
		double bias = 0.;
		auto n = cvs.size();

		// Reset vectors.
		std::fill(derivatives_.begin(), derivatives_.end(), 0);
		
		// Look up and apply grid bias. 
		if(grid_ != nullptr)
		{
			bool inbounds = true;
			std::vector<double> val(n, 0.);
			for(size_t i = 0; i < n; ++i)
			{
				val[i] = cvs[i]->GetValue();
				if(val[i] < grid_->GetLower(i) || val[i]  > grid_->GetUpper(i))
					inbounds = false;
			}

			if(inbounds)
			{
				auto frc = (*grid_)[val];
				for(size_t i = 0; i < n; ++i)
					derivatives_[i] = frc[i];
			}
			else
			{
				if(comm_.rank() == 0)
				{
					std::cerr << "Metadynamics: out of bounds ( ";
					for(auto& v : val)
						std::cerr << v << " "; 
					std::cerr << ")" << std::endl;
				}
			}
		}
		else
		{
			// Loop through hills and calculate the bias force.
			for(auto& hill : hills_)
			{		
				auto tbias = 1.;
				std::fill(tder_.begin(), tder_.end(), 1.0);
				std::fill(dx_.begin(), dx_.end(), 1.0);
				
				for(size_t i = 0; i < n; ++i)
				{
					dx_[i] = cvs[i]->GetDifference(hill.center[i]);
					tbias *= gaussian(dx_[i], hill.width[i]);
				}

				for(size_t i = 0; i < n; ++i)
					for(size_t j = 0; j < n; ++j)
					{
						if(j != i) 
							tder_[i] *= gaussian(dx_[j], hill.width[j]);
						else
							tder_[i] *= gaussianDerv(dx_[j], hill.width[j]);
					}

				bias += height_ * tbias;
				for(size_t i = 0; i < n; ++i)
					derivatives_[i] += height_*tder_[i];
			}
		}

		// Restraints.
		for(size_t i = 0; i < n; ++i)
		{
			auto cval = cvs[i]->GetValue();
			if(cval < lowerb_[i])
				derivatives_[i] += lowerk_[i]*(cval - lowerb_[i]);
			else if(cval > upperb_[i])
				derivatives_[i] += upperk_[i]*(cval - upperb_[i]);
		}
	}

}
