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

#ifndef _XML_CONFIG_FILE_H
#define _XML_CONFIG_FILE_H

#include <string>
#include <vector>
#include <iostream>
#include <map>

#include "XMLfile.h"
#include "streamline.h"

#include "linked_langevin_simulation.h"
#include "linked_moments_simulation.h"
#include "linked_pbm_simulation.h"
#include "stand_alone_langevin_simulation.h"
#include "stand_alone_moments_simulation.h"
#include "stand_alone_pbm_simulation.h"
#include "linked_moments_temp_simulation.h"
#include "linked_pbm_temp_simulation.h"

class XMLconfigfile : public XMLfile {

	/// Enumerator for Streamlines XML TAGS
	enum Config_Tags {
		NOT_DEF, 
			MESO_SIM, 
			CONDITIONS, 
				LINK_PATH, C_SPECIES,MOLAR_F_C, B_SPECIES, MOLAR_F_B, 
				PRESSURE, S_TEMP, E_TEMP, END_TIME, T_GRAD,
			DATA_PLOT,
				SAVE_STEPS, PRINT_STEPS, PARAMS, SAVE_PATH, SAVE_VTK_PATH, SAVE_VTK_STEPS,
			MODEL,
				MODEL_TYPE, VOLUME, MAX_PART, MIN_PART, MAX_TS, F_DIM, TIME_STEP
	};

	/// map for Strings and Tags
	std::map<std::string, Config_Tags> tag_map;

	/// Tags Map Initialization
	void map_Init();

	/// Returns data for a stand alone langevin simulation
	///	\param XMLElement* node: XML tree node
	raw_configuration_data get_stand_alone_langevin_data(XMLElement* node);

	/// Returns data for a stand alone pbm simulation
	///	\param XMLElement* node: XML tree node
	raw_configuration_data get_stand_alone_pbm_data(XMLElement* node);

	/// Returns data for a stand alone moments simulation
	///	\param XMLElement* node: XML tree node
	raw_configuration_data get_stand_alone_moments_data(XMLElement* node);

	/// Returns data for a linked langevin simulation
	///	\param XMLElement* node: XML tree node
	raw_configuration_data get_linked_langevin_data(XMLElement* node);

	/// Returns data for a linked moments simulation
	///	\param XMLElement* node: XML tree node
	raw_configuration_data get_linked_moments_data(XMLElement* node);

	/// Returns data for a linked pbm simulation
	///	\param XMLElement* node: XML tree node
	raw_configuration_data get_linked_pbm_data(XMLElement* node);

public:
	/// Constructor
	///	\param std::string path: xml file path
	XMLconfigfile(std::string _path);

	/// Runs the simulation specified in the XML file
	void run_simulation();


};

#endif // !_XML_CONFIG_FILE_H
