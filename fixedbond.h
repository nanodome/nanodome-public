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

#ifndef FIXEDBOND_H
#define FIXEDBOND_H

#include "pointbond.h"


template<typename P>
class FixedBond : public PointBond<P> {

    double bond_distance; ///< Distance between objects

	/// Link for constrains structures
	friend class Constrainer;

public:

    /// Constructor. Initialize the bond distance.
    FixedBond(std::shared_ptr<P> p0, std::shared_ptr<P> p1, double _d) : PointBond<P>(p0,p1), bond_distance(_d) {}

    /// Returns the bond distance which is fixed.
    double get_bond_distance() const { return bond_distance; }
};

#endif // FIXEDBOND_H
