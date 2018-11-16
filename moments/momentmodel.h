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

#ifndef MOMENTMODEL_H
#define MOMENTMODEL_H

#include "../cnt.h"
#include "../gasphase/gasphase.h"

/// Abstract class for the moments method models family.
class MomentModel {

protected:

    /// Condensing species, used for calculation of material properties:
    /// surface tension, saturation pressure, molecular mass, volume and surface.
    Species species;

public:

    /// Default constructor.
    MomentModel(Species _sp) : species(_sp) { }

    /// Nanoparticle density [#/m3]
    virtual double get_density() const = 0;

    /// Nanoparticle mean diameter [m]
    virtual double get_mean_diameter() const = 0;

    /// Total nanoparticles surface area [m2/m3]
    virtual double get_total_area() const = 0;

    /// Moment method timestep returning the condensing species consuption rate [#/m3 s]
    /// \param dt timestep size [s]
    /// \param T temperature [K]
    /// \param J nucleation rate [#/m3 s]
    /// \param j stable cluster size [#]
    /// \param S supersaturation ratio
    virtual double timestep(double dt, const GasPhase& gp, const NucleationTheory& nt) = 0;
};

#endif // MOMENTMODEL_H
