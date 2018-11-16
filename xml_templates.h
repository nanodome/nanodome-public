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

#ifndef XML_TEMPLATE_H
#define XML_TEMPLATE_H

#include <iostream>
#include <fstream>

void linked_moments_template(std::string _filename);

void linked_pbm_template(std::string _filename);

void linked_langevin_template(std::string _filename);

void stand_alone_moment_template(std::string _filename);

void stand_alone_pbm_template(std::string _filename);

void stand_alone_langevin_template(std::string _filename);




#endif
