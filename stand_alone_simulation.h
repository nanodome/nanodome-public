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

#ifndef STAND_ALONE_SIMULATION_H
#define STAND_ALONE_SIMULATION_H

#include "meso_simulation.h"
#include "gasphase/gasphasecv.h"

#include <list>

class StandAloneSimulation : public MesoSimulation {

protected:

	/// List of concentrating species
	std::list<std::string> c_species;

	/// list of concentrating species molar fraction (ordered with respect of c_species)
	std::list<double> c_species_m_fraction;

	/// List of bath species
	std::list<std::string> b_species;

	//list of bath species molar fraction(ordered with respect of b_species)
	std::list<double> b_species_m_fraction;

	/// Pressure [Pa]
	double pressure;

	/// Start Temperature [K]
	double start_T;

	/// End Temperature [K]
	double end_T;

	/// End Time [secs.]
	double end_Time;

	/// Temperature gradient [K/s]
	double temp_gradient;

public:

	/// Constructor
	///	\param: raw_configuration_data _xml_data raw data extracted from the XML configuration file
	StandAloneSimulation(raw_configuration_data _xml_data);

};


#endif
