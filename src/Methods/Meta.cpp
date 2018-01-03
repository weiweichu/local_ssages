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
#include "Meta.h"
#include "StandardMeta.h"
#include "WellTemperedMeta.h"
#include <math.h>
#include <iostream>
#include "Drivers/DriverException.h"
#include "CVs/CVManager.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"

using namespace Json;

namespace SSAGES
{
	// Post-integration hook.
	void Meta::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		auto cvs = cvmanager.GetCVs(cvmask_);
		// Add hills when needed.
		if(snapshot->GetIteration() % hillfreq_ == 0)
			AddHill(cvs, snapshot->GetIteration());

		// Always calculate the current bias.
		CalcBiasForce(cvs);
    if(grid_ != nullptr && potential_freq_ != 0 && snapshot->GetIteration() % potential_freq_ == 0)
    {
      PrintPotential(hills_.back(), snapshot->GetIteration());
    }
		// Take each CV and add its biased forces to the atoms
		// using the chain rule.
		auto& forces = snapshot->GetForces();
		auto& virial = snapshot->GetVirial();

		for(size_t i = 0; i < cvs.size(); ++i)
		{
			auto& grad = cvs[i]->GetGradient();
			auto& boxgrad = cvs[i]->GetBoxGradient();
			
			// Update the forces in snapshot by adding in the force bias from each
			// CV to each atom based on the gradient of the CV.
			for (size_t j = 0; j < forces.size(); ++j)
				for(size_t k = 0; k < 3; ++k)
					forces[j][k] -= derivatives_[i]*grad[j][k];
			
			virial += derivatives_[i]*boxgrad;
		}
	}

	// Post-simulation hook.
	void Meta::PostSimulation(Snapshot*, const CVManager&)
	{
	}

	// Load hills from file. 
	void Meta::LoadHills(const std::string& filename)
	{
		std::ifstream file(filename);
		std::string line; 

		auto dim = widths_.size();
		double iteration = 0, height=0;
		std::vector<double> width(dim, 0.), center(dim, 0);

		// First line is a comment. 
		std::getline(file, line);
		while(std::getline(file, line))
		{
			std::istringstream iss(line);
			iss >> iteration; // Not really using this. 
			
			// Load centers.
			for(size_t i = 0; i < dim; ++i)
				iss >> center[i];
			
			// Load widths.
			for(size_t i = 0; i < dim; ++i)
				iss >> width[i];
			
      iss >> height;

			hills_.emplace_back(center, width, height);
		}
	}


	// Writes hill to output file. This should only be called by the 
	// world master node. 
	void Meta::PrintHill(const Hill& hill, int iteration)
	{
		hillsout_.open("hills.out", std::fstream::app);
		
		hillsout_ << iteration + setoff_  << " ";
		hillsout_.precision(8);
		
		for(auto& cv : hill.center)
			hillsout_ << cv << " ";
		
		for(auto& w : hill.width)
			hillsout_ << w << " ";
    
		hillsout_ << hill.height << " ";
    hillsout_ << std::endl;
		hillsout_.close();
	}
  
  void Meta::PrintPotential(const Hill& hill, int iteration)
  {
   // std::cerr << "printpotential in meta " << iteration << "   " << hill.height << std::endl;
    potential_out_.open("potential.out", std::fstream::app);
    //potential_out_ << "timestep: " << iteration << "   " << "hill height: " << hill.height << std::endl;
    potential_out_ << iteration + setoff_ << "  " << hill.height << std::endl;
    auto it = grid_potential_->begin();
    for (; it != grid_potential_->end(); ++it)
    {
      auto& val = *it;
      potential_out_ << val << "   ";
    }
    potential_out_ << std::endl;
    potential_out_.close();
  }

  void Meta::ReadPotential(int timestep) {
    potential_initial_ = true;
    std::ifstream read_potential;
    std::string line, line1, line2, line3;
    int iteration = 0;
    if (timestep < 0) {
      read_potential.open("potential.out");
      while(std::getline(read_potential, line1)) {
        line3 = line2;
        line2 = line1;
      }
      read_potential.close();
      std::istringstream iss(line3);
      std::cout << "setoff:  " << setoff_ <<"   height: " << height_ << std::endl;
      iss >> setoff_;
      iss >> height_;
      std::cout << "setoff:  " << setoff_ <<"   height: " << height_ << std::endl;
      iss.clear();
      iss.str(line2);
      auto it = grid_potential_->begin();
      for(; it != grid_potential_->end(); ++it)
      {
        auto& val =  *it;
        iss >> val;
      }
    }
    else {
      read_potential.open("potential_pre.out");
      hillsout_.open("hills.out");
      hillsout_.close();
      while (std::getline(read_potential, line)) {
        std::istringstream iss(line);
        iss >> iteration;
        if (iteration == timestep) {
          iss >> height_;
          setoff_ = timestep;
          getline(read_potential, line);
          read_potential.close();
          std::istringstream iss(line);
          auto it = grid_potential_->begin();
          for(; it != grid_potential_->end(); ++it)
          {
            auto& val =  *it;
            iss >> val;
          }
          potential_out_.open("potential.out");
          potential_out_ << "step" << "   " << "height" << std::endl;
          potential_out_ << setoff_ << "   " << height_ << std::endl;
          it = grid_potential_->begin();
          for (; it != grid_potential_->end(); ++it)
          {
            auto& val = *it;
            potential_out_ << val << "   ";
          }
          potential_out_ << std::endl;
          break;
        }
      }
    }
  }
	Meta* Meta::Build(const Json::Value& json, 
		                  const MPI_Comm& world,
		                  const MPI_Comm& comm,
					      const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;
		
		reader.parse(JsonSchema::MetadynamicsMethod, schema);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		std::vector<double> widths;
		for(auto& s : json["widths"])
			widths.push_back(s.asDouble());
    
    auto height = json.get("height", 1.0).asDouble();
		std::vector<double> lowerb, upperb, lowerk, upperk;
		Grid<Vector>* grid = nullptr;
    Grid<double>* grid_potential = nullptr;
		if(json.isMember("grid")) 
    {
			grid = Grid<Vector>::BuildGrid(json.get("grid", Json::Value()));
      grid_potential = Grid<double>::BuildGrid(json.get("grid", Json::Value()));
    }
		else if(!json.isMember("lower_bounds") || !json.isMember("upper_bounds"))
			throw BuildException({
				"#/Method/Metadynamics: Both upper_bounds and lower_bounds "
				"must be defined if grid is not being used."});

		// Assume all vectors are the same size. 
		for(int i = 0; i < json["lower_bound_restraints"].size(); ++i)
		{
			lowerk.push_back(json["lower_bound_restraints"][i].asDouble());
			upperk.push_back(json["upper_bound_restraints"][i].asDouble());
			lowerb.push_back(json["lower_bounds"][i].asDouble());
			upperb.push_back(json["upper_bounds"][i].asDouble());
		}	
		auto hillfreq = json.get("hill_frequency", 1).asInt();
    auto potential_freq = json.get("potential_freq", 0).asInt();
		auto freq = json.get("frequency", 1).asInt();
    auto flavor = json.get("flavor", "none").asString();
    if (flavor == "WellTemperedMeta") {
      if(!json.isMember("grid")) {
        throw BuildException({
            "#/Method/Metadynamics: Grid "
            "must be defined if WellTemperedMeta flavor is being used."});
      }
      if(!json.isMember("delta_T")) {
        throw BuildException({
            "#/Method/Metadynamics: delta_T "
            "must be defined if WellTemperedMeta flavor is being used."});
      }
    //  grid_potential = Grid<double>::BuildGrid(json.get("grid", Json::Value()));
      auto delta_T = json.get("delta_T", 0.0).asDouble();
      auto* m = new WellTemperedMeta(
          world, comm, height, widths, 
          lowerb, upperb, lowerk, upperk,
          grid, grid_potential, hillfreq, potential_freq, freq, delta_T
          );

      if(json.isMember("load_hills")) {
      	m->LoadHills(json["load_hills"].asString());
      }
      if(json.isMember("read_potential")) {
        if(!json.isMember("grid")) {
          throw BuildException({
              "#/Method/Metadynamics: read_potential "
              "grid must be defined if read_potential is used"});
        }
        std::cerr << "start to read potential" << std::endl;
        auto timestep = json.get("read_potential", -1).asInt();
        //std::cerr << "potential_initial: " <<potential_initial_ << std::endl;
        m->ReadPotential(timestep);
        //std::cerr << "potential_initial: " <<potential_initial_ << std::endl;
      }
      return m;
    }
    
    else if (flavor == "Standard") {
      auto* m = new StandardMeta(
      world, comm, height, widths, 
      lowerb, upperb, lowerk,	upperk,
      grid, grid_potential, hillfreq, potential_freq, freq
      );

      if(json.isMember("load_hills"))
      	m->LoadHills(json["load_hills"].asString());
      if(json.isMember("read_potential")) {
        if(!json.isMember("grid")) {
          throw BuildException({
              "#/Method/Metadynamics: read_potential "
              "grid must be defined if read_potential is used"});
        }
        auto timestep = json.get("read_potential", 0).asInt();
        m->ReadPotential(timestep);
      }
      return m;
    }

    else {
      throw BuildException({"Unknown flavor for metadynamics. The options are \"WellTemperedMeta\" or \"Default\""});
    }
	}
}
