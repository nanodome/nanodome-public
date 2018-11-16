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

#ifndef CNT_H
#define CNT_H

#include "nucleationtheory.h"
#include "species.h"


/// Classical Nucleation Theory (CNT) implementation. This class is the implementation of the
/// CNT whose rates are calculated using the properties of the condensing species.
class ClassicalNucleationTheory : public NucleationTheory {

    Species species; ///< The nucleating species.

public:

    /// Default constructor.
    /// \param _species species type
    ClassicalNucleationTheory(Species _species);

    /// Primary particles formation rate [#/m3 s]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    double nucleation_rate(double T, double S) const;

    /// Smallest stable cluster particle number [#]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    double stable_cluster_size(double T, double S) const;

    /// Smallest stable cluster diameter [m]
    /// \param ns nucleating species concentration [#/m3]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    double stable_cluster_diameter(double T, double S) const;

    /// Surface condensation rate [#/m2 s]
    /// \param T temperature [K]
    /// \param S supersaturation ratio
    double condensation_rate(double T, double S) const;
};

#endif // CNT_H
