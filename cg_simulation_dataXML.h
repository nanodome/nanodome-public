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

#ifndef CGSIMULATION_DATAXML_H
#define CGSIMULATION_DATAXML_H

#include "XMLfile.h"
#include "aggregate/rattleaggregate.h"
#include "dynamicparticle.h"

#include <valarray>
#include <vector>
#include <memory>

class CGSimulationDataXML : public XMLfile {

	// RAW DATA extracted from the .vtp XML file

	/// raw data for the particles masses
	std::vector<double> raw_masses;

	/// raw data for diameters
	std::vector<double> raw_diameters;

	/// raw data for positions
	std::vector<double> positions;

	/// raw data for aggregates IDs
	std::vector<int> agg_IDs;

	/// raw data for particles IDs
	std::vector<int> part_IDs;

	/// raw data for connectivity
	std::vector<int> connections;

	/// raw data for edge types
	std::vector<int>edges_types;

	/// Number of aggregates in the simulation
	int aggregates_number;

	/// raw data parsing
	void get_raw_data();

public:

	/// Constructor
	///	\param std::string file_path path to the XML file
	CGSimulationDataXML(std::string file_path);

	/// Creates RATTLEAggregates from raw data
	std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>>
		create_aggregates();



};


#endif
