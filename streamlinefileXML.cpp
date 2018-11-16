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

#include "streamlinefileXML.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <iomanip>


std::vector<double> get_samples(const std::string& p_pcstStr, char delim);
std::list<Species> get_species(std::vector<std::string> s_species);
std::vector<std::string> tokenizer(const std::string& p_pcstStr, char delim);

std::string Vect_To_String(std::vector<double> vect);
std::string Species_To_String(std::list<Species> vect);


streamlinefileXML::streamlinefileXML(std::string _path, std::vector<streamline>& _streams) :
streams(_streams), XMLfile(_path)
{
	
	/// Get XML File
	std::cout<<"File path: "<<_path<<std::endl;
	path = _path;
	
	map_Init();

}

void streamlinefileXML::map_Init(){

    tag_map["GP"]               = GP;
    tag_map["T_START"]          = T_START;
    tag_map["T_END"]            = T_END;
    tag_map["ID"]               = ID;
    tag_map["STREAM"]           = STREAM;
    tag_map["N_TIME_SAMPLES"]   = N_TIME_SAMPLES;
    tag_map["TIME_SAMPLES"]     = TIME_SAMPLES;
    tag_map["T"]                = T;
    tag_map["P"]                = P;
    tag_map["N_SPECIES"]        = N_SPECIES;
    tag_map["SPECIES"]          = SPECIES;
    tag_map["MOLAR_C"]          = MOLAR_C;
	tag_map["MOLAR_F"]			= MOLAR_F;
	tag_map["X"]				= X;
	tag_map["Y"]				= Y;
	tag_map["Z"]				= Z;
}

void streamlinefileXML::read_Streamlines(){

	// Check file loading
	XMLError err = doc.LoadFile(path.c_str());

	if (err != XML_SUCCESS)
	{
		
		std::string error = "File not foud at path: " + path + "\n";
		error_entry_blocking(error);
	}
	
	// cursors to navigate the tree (root(e_root), leaves(e_curssor))
	XMLElement *e_root, *e_cursor;

	// Elaboration status variable
    bool status = false;

	// Get root of the XML File
	e_root = e_cursor = doc.FirstChildElement();


    status = check_Tag(e_root, "GP");
	Error_Check(status, e_cursor->Value());


	// Get number of streamlines
    e_cursor = e_cursor->FirstChildElement();

    status = check_Tag(e_cursor, "N_STREAM");
	Error_Check(status, e_cursor->Value());

    e_cursor->QueryIntText(&n_streams);

    // Get start time
    e_cursor = e_cursor->NextSiblingElement();

    status = check_Tag(e_cursor, "T_START");
	Error_Check(status, e_cursor->Value());

	e_cursor->QueryDoubleText(&samples_start);

	// Get end time
    e_cursor = e_cursor->NextSiblingElement();

    status = check_Tag(e_cursor, "T_END");
	Error_Check(status, e_cursor->Value());

	e_cursor->QueryDoubleText(&samples_end);

    // Checking Times
    /*if(samples_end <= samples_start){
		freopen(errFile.c_str(), "a", stderr);
        std::cerr<<"Streamlines Sampling start and end time Inconsistency"<<std::endl;
		fclose(stderr);
		exit(0);
    }*/

    // Get Streams
    e_cursor = e_root = e_root->FirstChildElement("STREAM");
	status = check_Tag(e_cursor, "STREAM");
	status = check_Tag(e_root, "STREAM");

	if (!status){
		std::string error = path + " Streamlines XML file Parsing Error entering STREAM sub-tree\n";
		error_entry_blocking(error);
	}

	// Parsing all STREAMS sub-tree
    int s = 0;
    while (s < n_streams && status) {

        status = parse_Stream(e_cursor, e_root);
		if (!status){
			std::string error = path + " Streamlines XML file Parsing Error reading STREAM sub-tree\n";
			error_entry_blocking(error);
		}
        e_cursor = e_root = e_root->NextSiblingElement();
        s++;
    }

    return;

}

bool streamlinefileXML::parse_Stream(XMLElement *cursor_ptr, XMLElement *root_ptr){

	bool stream_ok = true;

	// Stream data
	int s_ID = -1;
	int time_samples = -1;
	int n_species = -1;
	std::string raw_value;

	// Temporary streamlines data structures
	std::vector<double> s_Temp; std::vector<double> s_Press;
	std::vector<double> s_Mol; std::vector<double> s_Mol_frac;
	std::vector<double> s_X; std::vector<double> s_Y; std::vector<double> s_Z;
	std::vector<double> s_Time; std::list<Species> s_Species;

	// Parse subtree
	cursor_ptr = cursor_ptr->FirstChildElement();

	while (cursor_ptr != nullptr){

		// Getting the XML Tag for choosing the operation 
		std::string s_tag = cursor_ptr->Name();

		switch (tag_map[s_tag])
		{
		case ID: // get Stream ID
			
			cursor_ptr->QueryIntText(&s_ID);
			break;

		case N_TIME_SAMPLES: // get stream time samples cardinality
			cursor_ptr->QueryIntText(&time_samples);
			break;

		case N_SPECIES: // get stream species cardinality
			cursor_ptr->QueryIntText(&n_species);
			break;

		case TIME_SAMPLES: // get time samples
			raw_value = cursor_ptr->GetText();
			s_Time = get_samples(raw_value, ' ');
			break;

		case T: // get Temperature samples
			raw_value = cursor_ptr->GetText();
			s_Temp = get_samples(raw_value, ' ');
			break;

		case P: // get Pressure samples
			raw_value = cursor_ptr->GetText();
			s_Press = get_samples(raw_value, ' ');
			break;

		case SPECIES: // get Species
			raw_value = cursor_ptr->GetText();
			s_Species = get_species(tokenizer(raw_value, ' '));
			break;

		case MOLAR_C: // get species molar concentration
			raw_value = cursor_ptr->GetText();
			s_Mol = get_samples(raw_value, ' ');
			break;

		case MOLAR_F: // get species molar fraction
			raw_value = cursor_ptr->GetText();
			s_Mol_frac = get_samples(raw_value, ' ');
			break;

		case X:
			raw_value = cursor_ptr->GetText();
			s_X = get_samples(raw_value, ' ');
			break;
		case Y:
			raw_value = cursor_ptr->GetText();
			s_Y = get_samples(raw_value, ' ');
			break;
		case Z:
			raw_value = cursor_ptr->GetText();
			s_Z = get_samples(raw_value, ' ');
			break;

		default:
			std::string error = "Stream sub-tree Parsing error, XML tag <" + s_tag + "> not Recognized in STREAM subtree\n";
			error_entry_not_blocking(error);
			return stream_ok;

		}

		cursor_ptr = cursor_ptr->NextSiblingElement();
	}
	
	// Create Streamline
	streamline s_streamline(s_ID, s_Temp, s_Press, s_Time, s_Species, s_Mol, s_Mol_frac, s_X, s_Y, s_Z);

	streams.push_back(s_streamline);

	return stream_ok;
}

void streamlinefileXML::write_Streamlines(std::vector<streamline>& _streamlines,
									   double _start_t, double _end_t){

	// Create Root Node
	XMLNode *root_pt = doc.NewElement("GP");
	XMLElement *cursor_pt = root_pt->ToElement();

	
	doc.InsertEndChild(cursor_pt);

	// Create N_STREAMS leaf
	Create_XML_Node(cursor_pt, "N_STREAM", std::to_string(_streamlines.size()));

	// Create T_START leaf
	Create_XML_Node(cursor_pt, "T_START", std::to_string(_start_t));

	// Create T_END leaf
	Create_XML_Node(cursor_pt, "T_END", std::to_string(_end_t));

	// Create STREAM sub-trees
	for (auto it = _streamlines.begin(); it != _streamlines.end(); ++it){
		write_Stream((*it), cursor_pt);
	}
	
	// Save XML FIle in the specified path
	XMLError err = doc.SaveFile(path.c_str());

	if (err != XML_SUCCESS){
		std::string error = "Error in saving XML streamlines file in path: " + path + "\n";
		error_entry_not_blocking(error);
	}

	return;
}

void streamlinefileXML::write_Streamlines(std::string _path){

	// Create Root Node
	XMLNode *root_pt = doc.NewElement("GP");
	XMLElement *cursor_pt = root_pt->ToElement();


	doc.InsertEndChild(cursor_pt);

	// Create N_STREAMS leaf
	Create_XML_Node(cursor_pt, "N_STREAM", std::to_string(streams.size()));

	// Create T_START leaf
	Create_XML_Node(cursor_pt, "T_START", std::to_string(samples_start));

	// Create T_END leaf
	Create_XML_Node(cursor_pt, "T_END", std::to_string(samples_end));

	// Create STREAM sub-trees
	for (auto it = streams.begin(); it != streams.end(); ++it){
		write_Stream((*it), cursor_pt);
	}

	// Save XML FIle in the specified path
	XMLError err = doc.SaveFile(_path.c_str());

	if (err != XML_SUCCESS){
		std::string error = "Error in saving XML streamlines file in path: " + path + "\n";
		error_entry_not_blocking(error);
	}

	return;

}

void streamlinefileXML::write_Stream(streamline stream, XMLElement *r_sub_tree){

	XMLElement *cursor_ptr = doc.NewElement("STREAM");

	// Write Stream ID
	
	Create_XML_Node(cursor_ptr, "ID", std::to_string(stream.get_ID()));

	// Write Time samples
	Create_XML_Node(cursor_ptr, "N_TIME_SAMPLES", std::to_string((int)stream.get_Time().size()));

	// Write Time Samples
	Create_XML_Node(cursor_ptr, "TIME_SAMPLES", Vect_To_String(stream.get_Time()));

	// Write Temperatures
	Create_XML_Node(cursor_ptr, "T", Vect_To_String(stream.get_Temp()));

	// Write Presssures
	Create_XML_Node(cursor_ptr, "P", Vect_To_String(stream.get_Press()));

	// Write Species
	std::list<Species> s_list = stream.get_Species();

	// Write Number of species
	Create_XML_Node( cursor_ptr, "N_SPECIES", std::to_string(stream.get_Species().size() ) );

	// Write Species formula
	Create_XML_Node( cursor_ptr, "SPECIES", Species_To_String( stream.get_Species() ) );

	//Write Molar Concentrations
	Create_XML_Node(cursor_ptr, "MOLAR_C", Vect_To_String(stream.get_Molar()));

	//Write Molar Concentrations
	Create_XML_Node(cursor_ptr, "MOLAR_F", Vect_To_String(stream.get_Molar_Fraction()));

	// Write X coordinates
	Create_XML_Node(cursor_ptr, "X", Vect_To_String(stream.get_x()));

	// Write Y coordinates
	Create_XML_Node(cursor_ptr, "Y", Vect_To_String(stream.get_y()));

	// Write Z coordinates
	Create_XML_Node(cursor_ptr, "Z", Vect_To_String(stream.get_z()));

	// Add the STREAM Sub-tree to the GP root
	r_sub_tree->InsertEndChild(cursor_ptr);
}



void streamlinefileXML::printstreamlines(){

    for(auto it = streams.begin(); it != streams.end(); ++it)
        (*it).printstream();

}

// AUXILIARY FUNCTIONS

std::vector<double> get_samples(const std::string& p_pcstStr, char delim)
{
    std::vector<double> tokens;
    std::stringstream   mySstream(p_pcstStr);
    std::string temp;

    while (getline(mySstream, temp, delim)){
		if (!temp.empty())
			tokens.push_back(atof(temp.c_str()));
    }
    return tokens;
}

std::list<Species> get_species(std::vector<std::string> s_species){

    std::list<Species> s_list;
    for(int i = 0; i < s_species.size(); i++){

        Species s(s_species[i]);
        s_list.push_back(s);
    }

    return s_list;
}

std::vector<std::string> tokenizer(const std::string& p_pcstStr, char delim)
{
    std::vector<std::string> tokens;
    std::stringstream mySstream(p_pcstStr);
    std::string temp;

	while (getline(mySstream, temp, delim)){
		if (!temp.empty())
		tokens.push_back(temp);
	}

    return tokens;
}

std::string Vect_To_String(std::vector<double> vect){
	
	std::string ret_str;

	for (auto it : vect){
		std::stringstream ss;
		ss << std::setprecision(20) << it;
		std::string part;
		ss >> part;
		ret_str += part + " ";
	}

	return ret_str;
}


std::string Species_To_String(std::list<Species> vect){

	std::string ret_str;
	for (auto it : vect){
		std::string part = it.get_formula();
		ret_str += part + " ";
	}

	return ret_str;

}
