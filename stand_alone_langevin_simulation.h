/*
    NanoDome - H2020 European Project NanoDome, GA n.646121
    (www.nanodome.eu, https://github.com/nanodome/nanodome)
    e-mail: Emanuele Ghedini, emanuele.ghedini@unibo.it

    Copyright (C) 2018  Alma Mater Studiorum - Università di Bologna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef STAND_ALONE_LANGEVIN_SIMULATION
#define STAND_ALONE_LANGEVIN_SIMULATION

#include "stand_alone_simulation.h"
#include "dynamicparticle.h"
#include "aggregate/rattleaggregate.h"
#include "particlephase/constrainedlangevinparticlephase.h"
#include <experimental/filesystem>

class StandAloneLangevinSimulation : public StandAloneSimulation {

	/// Control Volume dimension [m3]
	double V;

	/// Max number of particles [#]
	int max_particles;

	/// Min number of particles [#]
	int min_particles;

	/// Path to save the VTK files
	std::string VTK_PATH;

	/// Frequency in saving VTK files
	int SAVE_VTK;

	/// Max value for the timestep (avoiding instability) [sec.]
	double max_dt;

public:
	/// Constructor
	///	\param: raw_configuration_data _xml_data XML file's raw data.
	StandAloneLangevinSimulation(raw_configuration_data _xml_data);

	/// runs the simulation
	void run_simulation();

};

#endif
