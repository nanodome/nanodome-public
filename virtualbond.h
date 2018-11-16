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

#ifndef VIRTUALBOND_H
#define VIRTUALBOND_H

#include "pointbond.h"


template<typename P>
class VirtualBond : public PointBond<P> {

    double cosalpha; ///< the cosine of the angle between particle p0 and p1

    std::shared_ptr<P> pv; ///< point on the angle vertex

    std::shared_ptr<PointBond<P>> b0; ///< the edge between pv and v0
    std::shared_ptr<PointBond<P>> b1; ///< the edge between pv and v1

public:

    /// default constructor
    /// \param p0 the first vertex
    /// \param p1 the second vertex
    /// \param _pv the corner vertx
    /// \param _b0 the edge between pv and p0
    /// \param _b1 the edge between pv and p1
    VirtualBond(std::shared_ptr<P> p0, std::shared_ptr<P> p1, std::shared_ptr<P> _pv, std::shared_ptr<PointBond<P>> _b0, std::shared_ptr<PointBond<P>> _b1);

    /// return the distance constraint calculated using cosalpha, e0 and e1 [m]
    double get_bond_distance() const;

    /// recalculate the angle depending on the actual b0 and b1 lengths
    void reset();
};

#include "virtualbond.cpp"

#endif // VIRTUALBOND_H
