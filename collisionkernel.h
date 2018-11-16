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

#ifndef COLLISIONKERNEL_H
#define COLLISIONKERNEL_H

#include "nanodome.h"
#include "collisionalobject.h"

#include <cmath>

/// Free molecular regime binary kernel for two collisional objects [m3/s]
/// Friedlander, "Smoke, Dust and Haze", p.192, eq 7.17
/// \param a0 the first colliding object
/// \param a1 the second colliding object
/// \param T system temperature [K]
double FreeMolecularKernel(CollisionalObject& a0, CollisionalObject& a1, double T) {

    double m0 = a0.get_mass();
    double m1 = a1.get_mass();
    double d0 = a0.get_collision_diameter();
    double d1 = a1.get_collision_diameter();

    return sqrt(0.5*M_PI*K_BOL*T)*sqrt(1/m0 + 1/m1) * square(d0 + d1);
}


/// Continuum regime binary kernel for two collisional objects [m3/s]
/// Friedlander, "Smoke, Dust and Haze", p.192, eq 7.16
/// \param a0 the first colliding object
/// \param a1 the second colliding object
/// \param T system temperature [K]
/// \param mu dynamic viscosity [Pa s]
double ContinuumKernel(CollisionalObject& a0, CollisionalObject& a1, double T, double mu) {

    double d0 = a0.get_collision_diameter();
    double d1 = a1.get_collision_diameter();

    return 2.0*K_BOL*T/(3.0*mu) * (1/d0 + 1/d1) * (d0 + d1);
}

#endif // COLLISIONKERNEL_H
