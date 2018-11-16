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

#include "XMLfile.h"

XMLfile::XMLfile(std::string _path) :path(_path) {
}

void XMLfile::Create_XML_Node(XMLElement *root_e, std::string TAG, std::string text){

	XMLElement *p_new_elem = doc.NewElement(TAG.c_str());
	p_new_elem->SetText(text.c_str());
	root_e->InsertEndChild(p_new_elem);

}

bool XMLfile::check_Tag(XMLElement *e_ptr, std::string TAG){

	bool check = false;

	if (e_ptr == nullptr || strcmp(e_ptr->Name(), TAG.c_str()))
	{
		std::string error = "Tag " + TAG + " NOT FOUND" + "\n";
		error_entry_blocking(error);
		return check;
	}
	else{
		check = true;
		std::cout << "Tag: " << TAG << " FOUND" << std::endl;
		return check;
	}

}

bool XMLfile::check_son(XMLElement *e_ptr, std::string TAG) {

	bool check = false;
	
	if (e_ptr->FirstChildElement(TAG.c_str()) == nullptr)
	{
		std::string error = "Tag " + TAG + " NOT FOUND in " + e_ptr->Name() + " TAG sons"+"\n";
		error_entry_blocking(error);
		return check;
	}
	else {
		check = true;
		//std::cout << "Tag: " << TAG << " FOUND" << std::endl;
		return check;
	}


}

void XMLfile::Error_Check(bool status, std::string TAG){

	if (!status){

		std::string error = "file: " + path + " Streamlines XML file Parsing Error at tag: " + TAG + "\n";
		error_entry_blocking(error);
	}
}
