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

#ifndef MOMENTMODELFRIEDLANDER_H
#define MOMENTMODELFRIEDLANDER_H

#include "momentmodel.h"


/// Implementation of the Friendlander model
/// S.K. Friedlander, "Smoke, Dust and Haze", Chapter 10
class MomentModelFriedlander : public MomentModel {

    double M0; ///< Nanoparticles number density [#/m3]
    double M1; ///< Nanoparticles total diameter [m/m3]
    double M2; ///< Total nanoparticles area [m2/m3]

	// test
	double M1_cond_test;

public:

    MomentModelFriedlander(Species _sp);

    /// Nanoparticle density [#/m3]
    double get_density() const { return M0; }

    /// Nanoparticle mean diameter [m]
    double get_mean_diameter() const { return (M0>0.0) ? M1/M0 : 0.0; }

    /// Total nanoparticles surface area [m2/m3]
    double get_total_area() const { return M2; }

	/// Get M1 condensation value
	double get_cond_term()const { return M1_cond_test; }

    /// Moment method timestep returning the condensing species consuption rate [#/m3 s]
    /// \param dt timestep size [s]
    /// \param T temperature [K]
    /// \param J nucleation rate [#/m3 s]
    /// \param j stable cluster size [#]
    /// \param S supersaturation ratio
    double timestep(double dt, const GasPhase& gp, const NucleationTheory& nt);
};


#endif // MOMENTMODELFRIEDLANDER_H
