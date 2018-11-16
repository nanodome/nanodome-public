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

#ifndef SPATIALAGGREGATE_H
#define SPATIALAGGREGATE_H

#include "fractalaggregate.h"

#include <valarray>

/// An aggregate type for which the geometry is described by the particle positions.
/// A restructuring SHAKE algorithm is implemented for rearranging the particles positions
/// after a change in the bonds constraints (e.g. for sintering).
///
/// Concepts:
/// P --> PointParticle

template<typename P>
class SpatialAggregate : public FractalAggregate<P> {

public:

    /// Constructor.
    /// \param p0 pointer to the newly created particle.
    SpatialAggregate(std::shared_ptr<P> p0) : FractalAggregate<P>(p0) { }

	/// Constructor from particles and bond lists
	SpatialAggregate(std::list<std::shared_ptr<P>> _particles,
		std::list<std::shared_ptr<ParticleBond<P>>> _bonds) :
		FractalAggregate<P>(_particles, _bonds) { }

	/// Copy Constructor
	SpatialAggregate(const SpatialAggregate& _aggregate) : 
		FractalAggregate<P>(_aggregate) {}

    /// Calculate the center of mass position [m]
    std::valarray<double> get_center_of_mass() const;

    /// Aggregate fractal dimension
    double get_fractal_dimension() const;

    /// Get the collision diameter [m]
    double get_collision_diameter() const;

    /// Get the diameter of the smallest sphere (centered in the center of mass)
    /// that encompasses completely the aggregate [m]
    double get_enclosing_sphere_diameter() const;

    /// Check collision between two aggregates
    /// !!! ALTERNATIVE: FRIEND FUNCTION
    /// \param a0 Other colliding aggregate
    /// \param c_dist maximum colliding distance
	///	\param std::shared_ptr<P>& p0 pointer to the colliding particles in the aggregate a0
	///	\param std::shared_ptr<P>& p1 pointer to the colliding particles in the aggregate a1
    bool check_collision(std::shared_ptr<SpatialAggregate<P>> a0, double c_dist,
                         std::shared_ptr<P>& p0, std::shared_ptr<P>& p1);

	/// Check collision between two aggregates (periodic conditions)
	/// \param a0 Other colliding aggregate
	/// \param c_dist maximum colliding distance
	///	\param std::valarray<double> shift coordinates shift in case of periodic conditions
	///	\param std::shared_ptr<P>& p0 pointer to the colliding particles in the aggregate a0
	///	\param std::shared_ptr<P>& p1 pointer to the colliding particles in the aggregate a1
	bool check_collision_periodic(std::shared_ptr<SpatialAggregate<P>> a0, double c_dist, std::valarray<double> shift,
								  std::shared_ptr<P>& p0, std::shared_ptr<P>& p1);

    /// Aggregate radius of gyration [m]
    double get_radius_of_gyration() const;

    /// SHAKE algorithm for bonds constraints implementation
    void shake();
};

#include "spatialaggregate.cpp"

#endif // SPATIALAGGREGATE_H
