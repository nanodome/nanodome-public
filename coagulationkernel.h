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

#ifndef COAGULATIONKERNEL_H
#define COAGULATIONKERNEL_H

#include "collisionalobject.h"
#include "gasphase/gasphase.h"

#include <list>
#include <memory>
#include <utility>


/// The basic interface for all kernels that work on a list of collisional objects.
/// It is meant to provide a coagulation rate for the ensemble and to return a pair
/// of randomly chosen colliding objects, according to kernel type.
///
/// Concept
/// T --> CollisionalObject
template<typename T>
class CoagulationKernel {

public:

    /// Get the coagulation rate for a given set of collisional objects [m3/s]
    /// \param objects the collisional objects set
    /// \param gp gas phase
    virtual double get_coag_rate(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp) = 0;

    /// Get the pair of coagulating particles for a given set of collisional objects. The
    /// function returns false if the coagulation is fictious.
    virtual std::pair<std::shared_ptr<T>,std::shared_ptr<T>> get_coag_particles(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp, double R_coag) = 0;
};

#endif // COAGULATIONKERNEL_H
