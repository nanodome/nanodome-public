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

#include "cg_simulation_dataXML.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <iomanip>

std::vector<double> get_values(const std::string& p_pcstStr, char delim);
std::vector<int> get_values_i(const std::string& p_pcstStr, char delim);

CGSimulationDataXML::CGSimulationDataXML(std::string file_path):aggregates_number(0), XMLfile(file_path) {

	// Check file loading
	XMLError err = doc.LoadFile(file_path.c_str());

	if (err != XML_SUCCESS)
	{

		std::string error = "File not foud at path: " + file_path + "\n";
		error_entry_blocking(error);
	}

	// get data from XML
	get_raw_data();
}

void CGSimulationDataXML::get_raw_data() {

	// cursors to navigate the tree (root(e_root), leaves(e_curssor))
	XMLElement *e_root, *e_cursor;

	// Elaboration status variable
	bool status = false;

	// Get root of the XML File
	e_root = e_cursor = doc.FirstChildElement();

	status = check_Tag(e_root, "VTKFile");
	Error_Check(status, e_cursor->Value());


	// Get PolyData TAG
	e_cursor = e_cursor->FirstChildElement();

	status = check_Tag(e_cursor, "PolyData");
	Error_Check(status, e_cursor->Value());

	// Get Piece TAG
	e_cursor = e_cursor->FirstChildElement();

	status = check_Tag(e_cursor, "Piece");
	Error_Check(status, e_cursor->Value());

	// Get number of particles
	std::string s_n_part = e_cursor->Attribute("NumberOfPoints");
	// Get number of edges
	std::string s_n_edges = e_cursor->Attribute("NumberOfLines");

	int particles_n = std::atoi(s_n_part.c_str());
	int edges_n = std::atoi(s_n_edges.c_str());

	// Get Points TAG
	e_cursor = e_cursor->FirstChildElement();
	status = check_Tag(e_cursor, "Points");
	Error_Check(status, e_cursor->Value());

	// Get Positions
	XMLElement *p_x = e_cursor->FirstChildElement();
	status = check_Tag(p_x, "DataArray");
	Error_Check(status, p_x->Value());

	std::string s_positions = p_x->GetText();

	positions = get_values(s_positions, ' ');

	// Get PointData
	e_cursor = e_cursor->NextSiblingElement();
	status = check_Tag(e_cursor, "PointData");
	Error_Check(status, e_cursor->Value());

	XMLElement *p_data = e_cursor->FirstChildElement();

	// Get Masses
	std::string s_masses = p_data->GetText();
	raw_masses = get_values(s_masses, ' ');

	// Get	particles diameters
	p_data = p_data->NextSiblingElement();
	std::string s_diameters = p_data->GetText();
	raw_diameters = get_values(s_diameters, ' ');

	// Skip	aggregates diameter
	p_data = p_data->NextSiblingElement();

	// Get aggregates IDs
	p_data = p_data->NextSiblingElement();
	std::string s_aggIDs = p_data->GetText();
	agg_IDs = get_values_i(s_aggIDs, ' ');

	// Get particles IDs
	p_data = p_data->NextSiblingElement();
	std::string s_partIDs = p_data->GetText();
	part_IDs = get_values_i(s_partIDs, ' ');

	// Get Edges Types
	e_cursor = e_cursor->NextSiblingElement();
	status = check_Tag(e_cursor, "CellData");
	Error_Check(status, e_cursor->Value());

	p_data = e_cursor->FirstChildElement();
	std::string s_types = p_data->GetText();
	edges_types = get_values_i(s_types, ' ');

	// Get Connectoivity data
	e_cursor = e_cursor->NextSiblingElement();
	status = check_Tag(e_cursor, "Lines");
	Error_Check(status, e_cursor->Value());

	p_data = e_cursor->FirstChildElement();
	std::string s_connectivity = p_data->GetText();
	connections = get_values_i(s_connectivity, ' ');
	
	
	return;


}

std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>>
CGSimulationDataXML::create_aggregates() {

	std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> agg_list;

	std::vector<std::shared_ptr<DynamicParticle>> particle_pool;
	std::vector<std::list<std::shared_ptr<ParticleBond<DynamicParticle>>>> bond_pool;

	int start_idx, end_idx, i;
	start_idx = end_idx = i = 0;

	std::vector<int>s_indices;
	std::vector<int>e_indices;
	int conn_cursor;
	int agg_ID_cursor;
	int part_cursor;
	int mass_cursor;
	int edge_t_cursor;
	int pos_cursor;
	conn_cursor = agg_ID_cursor = part_cursor = mass_cursor = edge_t_cursor = pos_cursor = 0;

	// Create particles pool
	while (start_idx < agg_IDs.size()) {
		
		// detect Aggregate
		int ID = agg_IDs[i];
		while (agg_IDs[i] == ID) {
			i++; end_idx++;
			if (i == agg_IDs.size())
				break;
		}
		// Update indices lists
		s_indices.push_back(start_idx);
		e_indices.push_back(end_idx);

		// Update aggregates number
		aggregates_number++;

		// Once got the intervall, create the aggregate
		for (int j = start_idx; j < end_idx; j++) {
			//get position
			double x1 = positions[pos_cursor]; pos_cursor++;
			double x2 = positions[pos_cursor]; pos_cursor++;
			double x3 = positions[pos_cursor]; pos_cursor++;
			std::valarray<double> x = { x1, x2, x3 };
			// get mass
			double m = raw_masses[mass_cursor]; mass_cursor++;

			// Species (MONOSPECIES)
			Species si("Si");
			// Get number of molecules composing the particles
			double n_molecules = m / si.get_mass();
			// Create the pointer to the particle
			auto p = std::make_shared<DynamicParticle>(n_molecules, si, x);
			// push it in the pool
			particle_pool.push_back(p);
		}
		start_idx = end_idx;
		
	}

	// Create bonds Pool
	agg_ID_cursor = 0;
	edge_t_cursor = 0;
	conn_cursor = 0;
	for (agg_ID_cursor = 0; agg_ID_cursor < aggregates_number; agg_ID_cursor++) {

		std::list<std::shared_ptr<ParticleBond<DynamicParticle>>> b_list;

		//for (int c_i = conn_cursor; c_i < connections.size(); c_i += 2) {

		int a = connections[conn_cursor]; conn_cursor++;
		int b = connections[conn_cursor]; conn_cursor++;
		while ((a >= s_indices[agg_ID_cursor] && a <= e_indices[agg_ID_cursor]) &&
			((b >= s_indices[agg_ID_cursor] && b <= e_indices[agg_ID_cursor]))) {

			if (edges_types[edge_t_cursor] != 1)
			{
				auto p = std::make_shared<ParticleBond<DynamicParticle>>(particle_pool[a], particle_pool[b]);
				b_list.push_back(p);
				if (conn_cursor >= connections.size())
					break;
				a = connections[conn_cursor]; conn_cursor++;
				b = connections[conn_cursor]; conn_cursor++;
				edge_t_cursor++;
			}
			else
			{
				if (conn_cursor >= connections.size())
					break;
				a = connections[conn_cursor]; conn_cursor++;
				b = connections[conn_cursor]; conn_cursor++;
				edge_t_cursor++;
			}
		}
		bond_pool.push_back(b_list);
		conn_cursor -= 2;
	}

	// Create Aggregates
	for (int i = 0; i < aggregates_number; i++) {
		std::list<std::shared_ptr<DynamicParticle>> p_list;
		// Create the list of particles composing the aggregate from the particles pool
		for (int j = s_indices[i]; j < e_indices[i]; j++) {
			p_list.push_back(particle_pool[j]);
		}
		// Create new aggregate
		auto agg_p = std::make_shared<RATTLEAggregate<DynamicParticle>>(p_list, bond_pool[i]);
		// Add it to the list
		agg_list.push_back(agg_p);
	}

	return agg_list;

}

std::vector<double> get_values(const std::string& p_pcstStr, char delim)
{
	std::vector<double> tokens;
	std::stringstream   mySstream(p_pcstStr);
	std::string temp;

	while (getline(mySstream, temp, delim)) {
		if (!temp.empty())
			tokens.push_back(atof(temp.c_str()));
	}
	return tokens;
}

std::vector<int> get_values_i(const std::string& p_pcstStr, char delim)
{
	std::vector<int> tokens;
	std::stringstream   mySstream(p_pcstStr);
	std::string temp;

	while (getline(mySstream, temp, delim)) {
		if (!temp.empty())
			tokens.push_back(atoi(temp.c_str()));
	}
	return tokens;
}
