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

#ifndef NUCLEATIONTHEORY_H
#define NUCLEATIONTHEORY_H

#include "nanodome.h"

/// Abstract class for nucleation theories. This class is an interface
/// for all nucleation theories that will be used in the gas and particle
/// phases.
class NucleationTheory {

public:

    /// Primary particles formation rate [#/m3 s]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    virtual double nucleation_rate(double T, double S) const = 0;

    /// Smallest stable cluster particle number [#]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    virtual double stable_cluster_size(double T, double S) const = 0;

    /// Smallest stable cluster diameter [m]
    /// Must return 1 if S<1
    /// \param ns nucleating species concentration [#/m3]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    virtual double stable_cluster_diameter(double T, double S) const = 0;

    /// Surface condensation rate [#/m2 s]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    virtual double condensation_rate(double T, double S) const = 0;
};

#endif // NUCLEATIONTHEORY_H
