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

#ifndef LINK_H
#define LINK_H

#include "tetra_vertex.h"
//#include "constants.h"

#include <utility>


class Link {

protected:
	/// edge's index in the adjacency matrix
	std::pair<int, int> e_index;

	/// tuple of pointers to the edge's extemes 
	std::pair<Vertex*, Vertex*> edge;

	/// Flag for indicating if is in the first tetrahedron
	bool base;

	/// Debug variable
	double stable_lenght;

public:

	/// Default Constructor
	Link();

	/// Constructor, edge index and couple of extremes
	Link(std::pair<Vertex*, Vertex*> _edge);

	/// Returns edge extremes
	std::pair<Vertex*, Vertex*> get_edge() const { return edge; }

	/// Returns edge indices
	std::pair<int, int> get_indices() const { return e_index; }

	/// Returns the first edge's extreme
	Vertex* get_first() const { return edge.first; }

	/// Returns the second edge's extreme
	Vertex* get_second() const { return edge.second; }

	/// Sets base member
	void set_base(bool _b) { base = _b; }

	/// Sets extremes
	void set_extremes(std::pair<Vertex*, Vertex*> p) { edge = p; }

	/// Returns base member value
	bool is_base() const { return base; }

	/// Returns the lenght of the edge
	virtual double get_dist() = 0;

	/// Shake constraint evaluation
	double rearrange();

	/// Debug function 
	double get_stable_dist() { return stable_lenght; }

	/// Debug function
	void set_stable_dist(double _dist) { stable_lenght = get_dist(); }
};



#endif // !EDGE_H
							
