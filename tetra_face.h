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

#ifndef FACE_H
#define FACE_H



#include "tetra_vertex.h"
#include "edge.h"
//#include "constants.h"
#include <tuple>





class Tetrahedron; // forward declaration


// Faces counter	
static int N_FACES = 0;


class Face
{
	/// face's index
	int f_index;

	/// list of vertices composing the face
	std::tuple<Vertex*, Vertex*, Vertex*> face;

	/// Plan identified by the face data
	/// equation: ax + by + cz - d = 0 parameters
	double a, b, c, d;

	/// centre of the face
	std::valarray<double> centre;

	/// Tetrahedra in which the face is included
	std::list<Tetrahedron*> tetra;

public:

	/// Blank Constructor
	Face();

	/// Constructor, passing face index and list of vertices composing the face
	/// \param _face: tuple with the pointers to the vertices of the face
	/// \param _edges: list of the edges connecting the face vertices
	Face(std::tuple<Vertex*, Vertex*, Vertex*> _face);

	/// Get distance between the centre of the face and a vertex
	/// \param _v: pointer to the vertex 
	double get_distance(Vertex* _v);

	/// Add tetrahedron in which the face is
	/// \param _t: pointer to the tetrahedron
	void add_tetrahedron(Tetrahedron* _t) { tetra.push_back(_t); }

	/// Returns the tuple of pointers to vertices composing the face
	std::tuple<Vertex*, Vertex*, Vertex*> get_vertices() const { return face; }

	/// Checks if the vertes is in the face
	/// \param _v:	Pointer to vertex to check
	bool is_in_face(Vertex* _v) const {
		if (std::get<0>(face)->get_local_id() == _v->get_local_id() ||
			std::get<1>(face)->get_local_id() == _v->get_local_id() ||
			std::get<2>(face)->get_local_id() == _v->get_local_id())
			return true;
		else
			return false;
	}

	/// Returns the other two vertices excluded the one in input
	std::list<Vertex*> get_opposite(Vertex* _v);

	/// Returns the Face GLOBAL ID
	int get_face_id() const { return f_index; }

};





#endif // !FACE_H
