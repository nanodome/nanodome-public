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

#ifndef STREAMLINEFILEJSON_H
#define STREAMLINEFILEJSON_H

/*IMPORTANT!!! when you change the platform and the compiler remember that the preprocessor definition in JSON.h 
               can give problems */


#include "JSONFile.h"
#include "streamline.h"



class streamlinefileJSON :public JSONfile{

	

public:

	/// Contructor
	/// \param std::string _file: path to the JSON file
	streamlinefileJSON(std::string _file);

	/// Parse the JSON file
	std::vector<streamline> parse();

	/// Write the JSON file starting from a set of streamlines
	void write_streamlines(std::vector<streamline>& _streams);

};



#endif // STREAMLINEFILEJSON_H
