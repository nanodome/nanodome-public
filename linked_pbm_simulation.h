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

#ifndef PBM_LINKED_SIMULATION
#define PBM_LINKED_SIMULATION

#include "linked_simulation.h"
#include "particlephase/pbmfractalparticlephase.h"
#include "aggregate/pbmaggregate.h"
#include "particle.h"

class PBMLinkedSimulation :public LinkedSimulation {

	/// Control Volume dimension [m3]
	double V;

	/// Max number of particles [#]
	int max_particles;

	/// Min number of particles [#]
	int min_particles;

	/// Fractal Dimension
	double fractal_dimension;

public:
	///Constructor
	///	\param: raw_configuration_data _xml_data XML file's raw data.
	PBMLinkedSimulation(raw_configuration_data _xml_data);

	/// runs the simulation
	void run_simulation();

};

#endif
