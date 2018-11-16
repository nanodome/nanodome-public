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

#ifndef DYNAMICPOINT_H
#define DYNAMICPOINT_H

#include "point.h"
#include "physicalobject.h"
#include "ndm_random.h"


class DynamicPoint : public Point , virtual public PhysicalObject {

    std::valarray<double> v; ///< Object velocity [m/s]
    std::valarray<double> f; ///< Object force [N]

public:

    /// Constructor.
    /// \param _x initial object position [m]
    /// \param _diameter object diameter [m]
    DynamicPoint(const std::valarray<double>& _x={0,0,0},
                 const std::valarray<double>& _v={0,0,0}) :
        Point(_x), v(_v), f(0.,3) { }

    /// Get velocity [m/s]
    std::valarray<double> get_velocity() const { return v; }

    /// Get force [N]
    std::valarray<double> get_force() const { return f; }

    /// Set a random maxwellian velocity
    /// \param temperature [K]
    void init_maxwellian_v(double T);

    /// Adjust force
    /// \param _f force correction [N]
    void add_force(const std::valarray<double>& _f) { f += _f; }

    /// Langevin-Verlet position step
    /// \param dt timestep [s]
    /// \param gamma friction coefficient [1/s]
    /// \param T temperature [K]
    void langevin_verlet_x(double dt, double gamma, double T);

    /// Langevin-Verlet velocity step
    /// \param dt timestep [s]
    void langevin_verlet_v(double dt);

    /// Box side reflection (TO IMPLEMENT PERIODIC BOUNDARIES ALSO)
    /// !!! POSSIBLE ALTERNATIVE: FRIEND FUNCTION
    /// \param side lenght of the volume box [m]
    void wall_reflection(double side);
};

#endif // DYNAMICOBJECT_H
