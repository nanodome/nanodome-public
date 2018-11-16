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

#include "tetra_vertex.h"
#include <climits>



Vertex::Vertex() :
	x(NULL),
	x_old(NULL),
	p(NULL),
	local_id(0), global_id(0),
	visited(false), hops(INT_MAX), base(false), added(false) { }

Vertex::Vertex(std::shared_ptr< DynamicParticle >_p, int local_index):
local_id(local_index), visited(false), hops(INT_MAX), 
base(false), added(false), p(_p) {

	x = &_p->x;
	x_old = &_p->x0;
	global_id = _p->get_id();

}


//Vertex::Vertex(std::valarray<double> _x, int _global_index, int _local_index) :
//	x(_x), x_old({ 0.0, 0.0, 0.0 }),
//	local_id(_local_index), global_id(_global_index),
//	visited(false), hops(INT_MAX), base(false), added(false) {}

//Vertex::Vertex(std::valarray<double> _x) :
//	x(_x), x_old({ 0.0, 0.0, 0.0 }),
//	local_id(0), global_id(-1),
//	visited(false), hops(INT_MAX), base(false), added(false) {}

