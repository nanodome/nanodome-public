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

#ifndef _XMLFILE_H
#define _XMLFILE_H

#include "i_o.h"
#include "tinyxml2/tinyxml2.h"

using namespace tinyxml2;

class XMLfile :public i_o{

protected:
	/// XML file Path
	std::string path;

	/// XML file (from tinyxml)
	XMLDocument doc;
	
	/// Create Node for the XML File
	/// \param XMLElement *cursor_ptr: pointer for expanding the <TAG> sub_tree
	/// \param std::string TAG: Name of the node
	/// \param string: Value of the node
	void Create_XML_Node(XMLElement *root_e, std::string TAG, std::string text);

	/// Check if the requested XML TAG is correct
	/// \param XMLElement *e_ptr: Pointer to the tinyxml ElementXML to check
	/// \param std::string TAG
	bool check_Tag(XMLElement *e_ptr, std::string TAG);

	/// Check if TAG is a son of e_ptr
	/// \param XMLElement *e_ptr: Pointer to the tinyxml ElementXML to check
	/// \param std::string TAG
	bool check_son(XMLElement *e_ptr, std::string TAG);

	/// Print out a detalled error and stops the execution
	/// bool status: status varible from other processes
	/// std::string TAG: XML tag giving problems
	void Error_Check(bool status, std::string TAG);

public:
	XMLfile(std::string _path);

};

#endif
