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

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "tetra_link.h"

/// Class for dummy bonds among particles to maintain the rigid body motion by means of a constrained algorithm
class Constraint : public Link
{
    /// Cosine of the angle subtended to the Contraint
    double cos_angle;

    // Subtended Vertex
    Vertex* angle_v;

    // Pointers to other edges composing the angle
    Link* p0;
    Link* p1;

public:

    /// Blank Constructor
    Constraint();

    /// Constructor
    /// \param _edge:	pair of pointers to the edge's extremes
    /// \param _p0:		pointer to the first side of the angle subtended by the edge
    /// \param _p1:		pointer to the second side of the angle subtended by the edge
    /// \param _angle:	pointer to the vertex of the angle subtended by the edge
    Constraint(std::pair<Vertex*, Vertex*> _edge, Link* _p0, Link* _p1, Vertex* _angle);

    /// Constructor
    //Constraint(std::pair<Vertex*, Vertex*> _edge);

    /// Returns Cosine of the angle
    double get_angle() const { return cos_angle; }

    /// Sets the angle subtebded by the constraint
    void set_angle(Vertex* _v, Link* _p0, Link* _p1);

    /// Returns the distance
    double get_dist();
};

#endif // CONSTRAINT_H
