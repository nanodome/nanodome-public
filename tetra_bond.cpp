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

#include "tetra_bond.h"
#include <cfloat>


PhysicalLink::PhysicalLink():
    lenght(DBL_MAX), b(NULL)
{
    e_index = std::pair<int, int>(-1, -1);
    edge = std::pair<Vertex*, Vertex*>(new Vertex(), new Vertex());
    base = false;
}

PhysicalLink::PhysicalLink(std::pair<Vertex*, Vertex*> _edge, std::shared_ptr<ParticleBond<DynamicParticle>> _b) :
    lenght(_b->get_bond_distance())
{
	b = _b;
    edge = _edge;
    e_index = std::pair<int, int>(edge.first->get_local_id(), edge.second->get_local_id());
    base = false;
}
