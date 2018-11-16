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

#include "XML_config_file.h"

#include <sstream>
#include <iostream>
#include <fstream>

std::list<std::string> get_tokens(const std::string& p_pcstStr, char delim)
{
	std::list<std::string> tokens;
	std::stringstream mySstream(p_pcstStr);
	std::string temp;

	while (getline(mySstream, temp, delim)) {
		if (!temp.empty())
			tokens.push_back(temp);
	}

	return tokens;
}

std::list<double> get_fractions(const std::string& p_pcstStr, char delim)
{
	std::list<double> tokens;
	std::stringstream mySstream(p_pcstStr);
	std::string temp;

	while (getline(mySstream, temp, delim)) {
		if (!temp.empty())
			tokens.push_back(atof(temp.c_str()));
	}

	return tokens;
}

XMLconfigfile::XMLconfigfile(std::string _path) : XMLfile(_path) {

	// file loading with checks
	XMLError err = doc.LoadFile(path.c_str());

	if (err != XML_SUCCESS)
	{

		std::string error = "File not found at path: " + path + "\n";
		error_entry_blocking(error);
	}

	// initialize tag map
	map_Init();

}

void XMLconfigfile::map_Init() {

	tag_map["MESO_SIM"]			= MESO_SIM;
	tag_map["CONDITIONS"]		= CONDITIONS;
	tag_map["LINK_PATH"]		= LINK_PATH;
	tag_map["C_SPECIES"]		= C_SPECIES;
	tag_map["MOLAR_F_C"]		= MOLAR_F_C;
	tag_map["B_SPECIES"]		= B_SPECIES;
	tag_map["MOLAR_F_B"]		= MOLAR_F_B;
	tag_map["PRESSURE"]			= PRESSURE;
	tag_map["S_TEMP"]			= S_TEMP;
	tag_map["E_TEMP"]			= E_TEMP;
	tag_map["END_TIME"]			= END_TIME;
	tag_map["T_GRAD"]			= T_GRAD;
	tag_map["DATA_PLOT"]		= DATA_PLOT;
	tag_map["SAVE_STEPS"]		= SAVE_STEPS;
	tag_map["PRINT_STEPS"]		= PRINT_STEPS;
	tag_map["PARAMS"]			= PARAMS;
	tag_map["SAVE_PATH"]		= SAVE_PATH;
	tag_map["SAVE_VTK_PATH"]	= SAVE_VTK_PATH;
	tag_map["SAVE_VTK_STEPS"]	= SAVE_VTK_STEPS;
	tag_map["MODEL"]			= MODEL;
	tag_map["MODEL_TYPE"]		= MODEL_TYPE;
	tag_map["VOLUME"]			= VOLUME;
	tag_map["MAX_PART"]			= MAX_PART;
	tag_map["MIN_PART"]			= MIN_PART;
	tag_map["MAX_TS"]			= MAX_TS;
	tag_map["F_DIM"]			= F_DIM;
	tag_map["TIME_STEP"]		= TIME_STEP;

}

void XMLconfigfile::run_simulation() {

	// elaboration status 
	bool status = false;
	// Linked Simulation
	bool linked = false;
	// linked Simulation only with temperature gradient
	bool temp_grad = false;
	// simulation types enumerator
	enum Simulation_Types {Not_def, Langevin, PBM, Moments };
	// Map for switch
	std::map<std::string, Simulation_Types> simulation_map;
	// init map
	simulation_map["Langevin"] = Langevin;
	simulation_map["PBM"] = PBM;
	simulation_map["Moments"] = Moments;

	// Analyze XML content
	// cursors to navigate the tree (root(e_root), leaves(e_curssor))
	XMLElement *e_root, *e_cursor;

	// Get root of the XML File
	e_root = e_cursor = doc.FirstChildElement();

	// Checks if is a Nanodome COnfiguration file
	status = check_Tag(e_root, "MESO_SIM");
	Error_Check(status, e_cursor->Value());

	// Check if it is a linked simulation (checks if the LINK_PATH TAG exists)
	if (e_cursor->FirstChildElement("CONDITIONS")->FirstChildElement("LINK_PATH") != nullptr)
		linked = true;
	
	// Check if it is a linked simulatiion only with temperature gradient
	if (linked) {
		std::string t_grad = e_cursor->FirstChildElement("CONDITIONS")->FirstChildElement("LINK_PATH")->Attribute("temp_grad");
		if (t_grad == "yes")
			temp_grad = true;
	}

	//Get the simulation type
	raw_configuration_data sim_data;
	sim_data.model_type = e_cursor->FirstChildElement("MODEL")->FirstChildElement("MODEL_TYPE")->GetText();

	// Start proper simulation

	if (linked) {
		switch(simulation_map[sim_data.model_type]){
		case Langevin:
		{
			std::cout << "Langevin, Linked Meso simulation" << std::endl;
			sim_data = get_linked_langevin_data(e_cursor);
			LinkedLangevinSimulation sim(sim_data);
			sim.run_simulation();
		}
			break;
		case PBM:
		{
			std::cout << "PBM, Linked Meso simulation" << std::endl;
			sim_data = get_linked_pbm_data(e_cursor);
			if (temp_grad) {
				PBMLinkedTempSimulation sim(sim_data);
				sim.run_simulation();
			}
			else {
				PBMLinkedSimulation sim(sim_data);
				sim.run_simulation();
			}
		}
			break;
		case Moments: {
			std::cout << "Moments, Linked Meso simulation" << std::endl;
			sim_data = get_linked_moments_data(e_cursor);
			if (temp_grad) {
				LinkedMomentsTempSimulation sim(sim_data);
				sim.run_simulation();
			}
			else {
				LinkedMomentsSimulation sim(sim_data);
				sim.run_simulation();
			}
			
		}
			break;

		default:
			std::cout << "Simulation not recognized!! QUITTING" << std::endl;
			break;
		}
	}
	else {
		switch (simulation_map[sim_data.model_type]) {
		case Langevin: {
			std::cout << "Langevin, Stand Alone Meso simulation" << std::endl;
			sim_data = get_stand_alone_langevin_data(e_cursor);
			StandAloneLangevinSimulation sim(sim_data);
			sim.run_simulation();
		}
			break;
		case PBM:
		{
			std::cout << "PBM, Stand Alone Meso simulation" << std::endl;
			sim_data = get_stand_alone_pbm_data(e_cursor);
			StandAlonePBMSimulation sim(sim_data);
			sim.run_simulation();
		}
			break;
		case Moments: 
		{
			std::cout << "Moments, Stand Alone Meso simulation" << std::endl;
			sim_data = get_stand_alone_moments_data(e_cursor);
			StandAloneMomentsSimulation sim(sim_data);
			sim.run_simulation();
		}
			break;
		default:
			std::cout << "Simulation not recognized!! QUITTING" << std::endl;
			break;
		}
		
	}
}

raw_configuration_data XMLconfigfile::get_stand_alone_langevin_data(XMLElement* node) {

	raw_configuration_data ret;

	// Get simulation boundary condition
	// Condensating species
	if (check_son(node->FirstChildElement("CONDITIONS"), "C_SPECIES")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("C_SPECIES")->GetText();
		ret.c_species = get_tokens(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG C_SPECIES");
	}
	// Carrier gas species
	if (check_son(node->FirstChildElement("CONDITIONS"), "B_SPECIES")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("B_SPECIES")->GetText();
		ret.b_species = get_tokens(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG B_SPECIES");
	}

	// Condensating species molar fraction
	if (check_son(node->FirstChildElement("CONDITIONS"), "MOLAR_F_C")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("MOLAR_F_C")->GetText();
		ret.c_s_molar_fraction = get_fractions(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG MOLAR_F_C");
	}

	// carrier gas species molar fraction
	if (check_son(node->FirstChildElement("CONDITIONS"), "MOLAR_F_B")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("MOLAR_F_B")->GetText();
		ret.b_s_molar_fraction = get_fractions(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG MOLAR_F_B");
	}

	// Pressure
	if (check_son(node->FirstChildElement("CONDITIONS"), "PRESSURE")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("PRESSURE")->GetText();
		ret.pressure = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG PRESSURE");
	}

	// Start Temperature
	if (check_son(node->FirstChildElement("CONDITIONS"), "S_TEMP")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("S_TEMP")->GetText();
		ret.start_temp = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG S_TEMP");
	}
	// End Temperature
	if (check_son(node->FirstChildElement("CONDITIONS"), "E_TEMP")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("E_TEMP")->GetText();
		ret.end_temp = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG E_TEMP");
	}
	// End Time
	if (check_son(node->FirstChildElement("CONDITIONS"), "END_TIME")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("END_TIME")->GetText();
		ret.end_time = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG END_TIME");
	}
	// Temperature Gradient
	if (check_son(node->FirstChildElement("CONDITIONS"), "T_GRAD")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("T_GRAD")->GetText();
		ret.t_grad = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG T_GRAD");
	}
	// Get Plot data
	// Frequency in saving on file data every SAVE_STEP steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_STEPS")->GetText();
		ret.SAVE_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_STEPS");
	}
	// Frequency in printing on screen file data every PRINT_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "PRINT_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("PRINT_STEPS")->GetText();
		ret.PRINT_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG PRINT_STEPS");
	}
	// Frequency in saving VTK visualization file each SAVE_VTK_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_VTK_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_VTK_STEPS")->GetText();
		ret.SAVE_VTK_STEPS = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_VTK_STEPS");
	}
	// Path to save VTK data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_VTK_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.vtk_save = sub_tree->FirstChildElement("SAVE_VTK_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_VTK_PATH");
	}
	// path for saving data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.save_path = sub_tree->FirstChildElement("SAVE_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_PATH");
	}
	// get Simulation data
	// Volume
	if (check_son(node->FirstChildElement("MODEL"), "VOLUME")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("VOLUME")->GetText();
		ret.volume = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG VOLUME");
	}
	// Max number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MAX_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MAX_PART")->GetText();
		ret.max_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MAX_PART");
	}
	// Min number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MIN_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MIN_PART")->GetText();
		ret.min_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MIN_PART");
	}
	// max time step 
	if (check_son(node->FirstChildElement("MODEL"), "MAX_TS")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MAX_TS")->GetText();
		ret.max_time_step = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MAX_TS");
	}
	
	
	return ret;
}
raw_configuration_data XMLconfigfile::get_stand_alone_pbm_data(XMLElement* node) {

	raw_configuration_data ret;

	// Condensating species
	if (check_son(node->FirstChildElement("CONDITIONS"), "C_SPECIES")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("C_SPECIES")->GetText();
		ret.c_species = get_tokens(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG C_SPECIES");
	}
	
	// Carrier gas species
	if (check_son(node->FirstChildElement("CONDITIONS"), "B_SPECIES")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("B_SPECIES")->GetText();
		ret.b_species = get_tokens(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG B_SPECIES");
	}

	// Condensating species molar fraction
	if (check_son(node->FirstChildElement("CONDITIONS"), "MOLAR_F_C")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("MOLAR_F_C")->GetText();
		ret.c_s_molar_fraction = get_fractions(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG MOLAR_F_C");
	}

	// carrier gas species molar fraction
	if (check_son(node->FirstChildElement("CONDITIONS"), "MOLAR_F_B")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("MOLAR_F_B")->GetText();
		ret.b_s_molar_fraction = get_fractions(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG MOLAR_F_B");
	}

	// Pressure
	if (check_son(node->FirstChildElement("CONDITIONS"), "PRESSURE")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("PRESSURE")->GetText();
		ret.pressure = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG PRESSURE");
	}

	// Start Temperature
	if (check_son(node->FirstChildElement("CONDITIONS"), "S_TEMP")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("S_TEMP")->GetText();
		ret.start_temp = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG S_TEMP");
	}

	// End Temperature
	if (check_son(node->FirstChildElement("CONDITIONS"), "E_TEMP")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("E_TEMP")->GetText();
		ret.end_temp = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG E_TEMP");
	}

	// End Time
	if (check_son(node->FirstChildElement("CONDITIONS"), "END_TIME")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("END_TIME")->GetText();
		ret.end_time = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG END_TIME");
	}

	// Temperature Gradient
	if (check_son(node->FirstChildElement("CONDITIONS"), "T_GRAD")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("T_GRAD")->GetText();
		ret.t_grad = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG T_GRAD");
	}

	// Get Plot data
	// Frequency in saving on file data every SAVE_STEP steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_STEPS")->GetText();
		ret.SAVE_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_STEPS");
	}

	// Frequency in printing on screen file data every PRINT_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "PRINT_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("PRINT_STEPS")->GetText();
		ret.PRINT_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG PRINT_STEPS");
	}

	// path for saving data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.save_path = sub_tree->FirstChildElement("SAVE_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_PATH");
	}

	// get Simulation data
	// Volume
	if (check_son(node->FirstChildElement("MODEL"), "VOLUME")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("VOLUME")->GetText();
		ret.volume = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG VOLUME");
	}

	// Max number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MAX_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MAX_PART")->GetText();
		ret.max_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MAX_PART");
	}

	// Min number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MIN_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MIN_PART")->GetText();
		ret.min_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MIN_PART");
	}

	// Fractal Dimension
	if (check_son(node->FirstChildElement("MODEL"), "F_DIM")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("F_DIM")->GetText();
		ret.frac_dim = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG F_DIM");
	}


	return ret;
}
raw_configuration_data XMLconfigfile::get_stand_alone_moments_data(XMLElement* node) {
	
	raw_configuration_data ret;
	// Condensating species
	if (check_son(node->FirstChildElement("CONDITIONS"), "C_SPECIES")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("C_SPECIES")->GetText();
		ret.c_species = get_tokens(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG C_SPECIES");
	}

	// Carrier gas species
	if (check_son(node->FirstChildElement("CONDITIONS"), "B_SPECIES")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("B_SPECIES")->GetText();
		ret.b_species = get_tokens(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG B_SPECIES");
	}

	// Condensating species molar fraction
	if (check_son(node->FirstChildElement("CONDITIONS"), "MOLAR_F_C")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("MOLAR_F_C")->GetText();
		ret.c_s_molar_fraction = get_fractions(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG MOLAR_F_C");
	}

	// carrier gas species molar fraction
	if (check_son(node->FirstChildElement("CONDITIONS"), "MOLAR_F_B")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("MOLAR_F_B")->GetText();
		ret.b_s_molar_fraction = get_fractions(s_data, ' ');
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG MOLAR_F_B");
	}

	// Pressure
	if (check_son(node->FirstChildElement("CONDITIONS"), "PRESSURE")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("PRESSURE")->GetText();
		ret.pressure = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG PRESSURE");
	}

	// Start Temperature
	if (check_son(node->FirstChildElement("CONDITIONS"), "S_TEMP")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("S_TEMP")->GetText();
		ret.start_temp = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG S_TEMP");
	}

	// End Temperature
	if (check_son(node->FirstChildElement("CONDITIONS"), "E_TEMP")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("E_TEMP")->GetText();
		ret.end_temp = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG E_TEMP");
	}

	// End Time
	if (check_son(node->FirstChildElement("CONDITIONS"), "END_TIME")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("END_TIME")->GetText();
		ret.end_time = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG END_TIME");
	}

	// Temperature Gradient
	if (check_son(node->FirstChildElement("CONDITIONS"), "T_GRAD")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		std::string s_data = sub_tree->FirstChildElement("T_GRAD")->GetText();
		ret.t_grad = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG T_GRAD");
	}

	// Get Plot data
	// Frequency in saving on file data every SAVE_STEP steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_STEPS")->GetText();
		ret.SAVE_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_STEPS");
	}

	// Frequency in printing on screen file data every PRINT_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "PRINT_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("PRINT_STEPS")->GetText();
		ret.PRINT_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG PRINT_STEPS");
	}

	// path for saving data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.save_path = sub_tree->FirstChildElement("SAVE_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_PATH");
	}
	// get timestep
	if (check_son(node->FirstChildElement("MODEL"), "TIME_STEP")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("TIME_STEP")->GetText();
		ret.dt = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG TIME_STEP");
	}
	// Get Model type
	if (check_son(node->FirstChildElement("MODEL"), "MODEL_TYPE")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MODEL_TYPE")->GetText();
		ret.model_type = s_data.c_str();
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MODEL_TYPE");
	}

	return ret;
}
raw_configuration_data XMLconfigfile::get_linked_langevin_data(XMLElement* node) {

	raw_configuration_data ret;

	// get streamlines file's path
	if (check_son(node->FirstChildElement("CONDITIONS"), "LINK_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		ret.streams_path = sub_tree->FirstChildElement("LINK_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG LINK_PATH");
	}
	// Get Plot data
	// Frequency in saving on file data every SAVE_STEP steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_STEPS")->GetText();
		ret.SAVE_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_STEPS");
	}

	// Frequency in printing on screen file data every PRINT_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "PRINT_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("PRINT_STEPS")->GetText();
		ret.PRINT_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG PRINT_STEPS");
	}

	// Frequency in saving VTK visualization file each SAVE_VTK_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_VTK_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_VTK_STEPS")->GetText();
		ret.SAVE_VTK_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_VTK_STEPS");
	}

	// Path to save VTK data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_VTK_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.vtk_save = sub_tree->FirstChildElement("SAVE_VTK_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_VTK_PATH");
	}

	// path for saving data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.save_path = sub_tree->FirstChildElement("SAVE_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_PATH");
	}

	// get Simulation data
	// Volume
	if (check_son(node->FirstChildElement("MODEL"), "VOLUME")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("VOLUME")->GetText();
		ret.volume = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG VOLUME");
	}

	// Max number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MAX_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MAX_PART")->GetText();
		ret.max_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MAX_PART");
	}

	// Min number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MIN_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MIN_PART")->GetText();
		ret.min_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MIN_PART");
	}

	// max time step 
	if (check_son(node->FirstChildElement("MODEL"), "MAX_TS")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MAX_TS")->GetText();
		ret.max_time_step = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MAX_TS");
	}

	return ret;
}
raw_configuration_data XMLconfigfile::get_linked_moments_data(XMLElement* node) {

	raw_configuration_data ret;

	// get streamlines file's path
	if (check_son(node->FirstChildElement("CONDITIONS"), "LINK_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		ret.streams_path = sub_tree->FirstChildElement("LINK_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG LINK_PATH");
	}

	// Get Plot data
	// Frequency in saving on file data every SAVE_STEP steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_STEPS")->GetText();
		ret.SAVE_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_STEPS");
	}

	// Frequency in printing on screen file data every PRINT_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "PRINT_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("PRINT_STEPS")->GetText();
		ret.PRINT_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG PRINT_STEPS");
	}

	// path for saving data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.save_path = sub_tree->FirstChildElement("SAVE_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_PATH");
	}

	// get timestep
	if (check_son(node->FirstChildElement("MODEL"), "TIME_STEP")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("TIME_STEP")->GetText();
		ret.dt = atof(s_data.c_str());
	}
	
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG TIME_STEP");
	}

	// Get Model type
	if (check_son(node->FirstChildElement("MODEL"), "MODEL_TYPE")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MODEL_TYPE")->GetText();
		ret.model_type = s_data.c_str();
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MODEL_TYPE");
	}

	return ret;
}
raw_configuration_data XMLconfigfile::get_linked_pbm_data(XMLElement* node) {

	raw_configuration_data ret;

	// get streamlines file's path
	if (check_son(node->FirstChildElement("CONDITIONS"), "LINK_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("CONDITIONS");
		ret.streams_path = sub_tree->FirstChildElement("LINK_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In CONDITIONS SUBTREE MISSING TAG LINK_PATH");
	}
	// Get Plot data
	// Frequency in saving on file data every SAVE_STEP steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("SAVE_STEPS")->GetText();
		ret.SAVE_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_STEPS");
	}
	// Frequency in printing on screen file data every PRINT_STEPS steps
	if (check_son(node->FirstChildElement("DATA_PLOT"), "PRINT_STEPS")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		std::string s_data = sub_tree->FirstChildElement("PRINT_STEPS")->GetText();
		ret.PRINT_STEPS = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG PRINT_STEPS");
	}
	
	// path for saving data
	if (check_son(node->FirstChildElement("DATA_PLOT"), "SAVE_PATH")) {
		XMLElement* sub_tree = node->FirstChildElement("DATA_PLOT");
		ret.save_path = sub_tree->FirstChildElement("SAVE_PATH")->GetText();
	}
	else {
		i_o error;
		error.error_entry_blocking("In DATA_PLOT SUBTREE MISSING TAG SAVE_PATH");
	}
	// get Simulation data
	// Volume
	if (check_son(node->FirstChildElement("MODEL"), "VOLUME")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("VOLUME")->GetText();
		ret.volume = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG VOLUME");
	}
	// Max number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MAX_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MAX_PART")->GetText();
		ret.max_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MAX_PART");
	}
	// Min number of particles
	if (check_son(node->FirstChildElement("MODEL"), "MIN_PART")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("MIN_PART")->GetText();
		ret.min_part = atoi(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG MIN_PART");
	}
	
	// Fractal Dimension
	if (check_son(node->FirstChildElement("MODEL"), "F_DIM")) {
		XMLElement* sub_tree = node->FirstChildElement("MODEL");
		std::string s_data = sub_tree->FirstChildElement("F_DIM")->GetText();
		ret.frac_dim = atof(s_data.c_str());
	}
	else {
		i_o error;
		error.error_entry_blocking("In MODEL SUBTREE MISSING TAG F_DIM");
	}

	return ret;
}
