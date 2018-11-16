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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "mesoobject.h"
#include "objectcounter.h"
#include "species.h"

#include <memory>


/// Particle implements the data structure needed to describe a primary particle
/// and the basic functions to change its composition and to describe it.
/// A particle is a collection of molecules of different species. The
/// particle shape is assumed to be spherical.
class Particle : public MesoObject, public ObjectCounter<Particle> {

    double n;  ///< number of monomers
    Species s; ///< species composing the particle

public:

    /// Constructor.
    /// \param _n number of molecules in the particle [#]
    /// \param _s species type
    Particle(double _n, Species _s) : n(_n), s(_s) { }

	
    /// Copy constructor
    Particle(const Particle& p1);

    /// Get particle molecules number [#].
    double get_n() const { return n; }

    /// Get particle species.
    const Species& get_species() const { return s; }

    /// Add molecules.
    /// \param dn amount of molecules to add [#]
    void add_molecules(double dn) { n += dn; }

    /// Get particle mass from the molecules composition [kg]
    double get_mass() const;

    /// Get particle volume using the molecules volume [m3]
    double get_volume() const;

    /// Get particle diameter [m]
    double get_diameter() const;

    /// Get particle surface (spherical assumption) [m2]
    double get_surface() const;
};

#endif // PARTICLE_H
