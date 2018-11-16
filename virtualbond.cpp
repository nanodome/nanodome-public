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

#include "virtualbond.h"


template<typename P>
VirtualBond<P>::VirtualBond(std::shared_ptr<P> p0, std::shared_ptr<P> p1, std::shared_ptr<P> _pv, std::shared_ptr<PointBond<P>> _b0, std::shared_ptr<PointBond<P>> _b1)
    : PointBond<P>(p0,p1), pv(_pv), b0(_b0), b1(_b1)
{
    reset();
}


template<typename P>
double VirtualBond<P>::get_bond_distance() const {

    // get the e0 and e1 distance constraints
    double d0 = b0->get_bond_distance();
    double d1 = b1->get_bond_distance();

    // cosine law
    return std::sqrt(d0*d0 + d1*d1 - 2*d0*d1*cosalpha);
}


template<typename P>
void VirtualBond<P>::reset() {

    // calculate the angle
    std::valarray<double> rv0 = this->get_v0()->get_x() - pv->get_x();
    std::valarray<double> rv1 = this->get_v1()->get_x() - pv->get_x();

    cosalpha = (rv0*rv1).sum()/(sqrt((rv0*rv0).sum()) * sqrt((rv1*rv1).sum()));
}
