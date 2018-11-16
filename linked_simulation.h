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

#ifndef LINKED_SIMULATION_H
#define LINKED_SIMULATION_H

#include "meso_simulation.h"
#include "streamlinefileXML.h"
#include "utilities.h"
#include "gasphase/gasphaselink.h"

#include <vector>

class LinkedSimulation :public MesoSimulation {

protected:

	/// path to store the streamlines data
	std::string SAVE_STREAMS;

	/// path for retreiving the XML file with the streamlines
	std::string STREAMS_PATH;

	/// List of streamlines in the linked simulation
	std::vector<streamline> l_streams;

public:
	/// Constructor
	///	\param: raw_configuration_data _xml_data
	LinkedSimulation(raw_configuration_data _xml_data);


};



#endif
