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

#ifndef MOMENTMODELPRATSINIS_H
#define MOMENTMODELPRATSINIS_H

#include "momentmodel.h"


/// Implementation of the Pratsinis 1988 model.
/// Simultaneous Nucleation, Condensation and Coagulation in Aerosol Reactors,
/// S.E. Pratsinis, Journal of Colloid and Interface Science, Vol. 124, No. 2, August 1988
class MomentModelPratsinis :public MomentModel {

    double M0; ///< Nanoparticles number density [#/m3]
    double M1; ///< Nanoparticles total volume [m3/m3]
    double M2; ///< Nanoparticles ... [m6/m3]

    // auxiliary functions defined and used by the model
    double M_k(double k, double sg, double vg);

    double zeta0(double sg);
    double zeta2(double sg);

    double csi1(double T);
    double csi2(double T);

    double ln2_standard_dev();

	/// TEST
	double M1_cond;
	

public:

    /// Standard constructor.
    MomentModelPratsinis(Species _sp);

	/// Constructor with initialization for the moments
	MomentModelPratsinis(Species _sp, double _M0, double _M1, double _M2);

    /// Nanoparticle density [#/m3]
    double get_density() const { return M0; }

    /// Nanoparticle mean diameter [m]
    double get_mean_diameter() const { return (M0>0.0) ? pow(6*M1/(M_PI*M0),1./3.) : 0.0; }

    /// Nanoparticle volume density [m3/m3]
    double get_total_volume() const { return M1; }

    /// ... [m6/m3]
    double get_total_area() const { return M2; }

    /// Timestep calculation
    /// \param dt timestep size [s]
    /// \param T temperature [K]
    /// \param J nucleation rate [#/m3 s]
    /// \param j stable cluster size [#]
    /// \param S supersaturation ratio
    double timestep(double dt,  const GasPhase& gp, const NucleationTheory& nt);

	/// Timestep calculation
	/// \param dt timestep size [s]
	/// \param T temperature [K]
	/// \param J nucleation rate [#/m3 s]
	/// \param j stable cluster size [#]
	/// \param S supersaturation ratio
	double timestep(double dt, const GasPhase& gp, const NucleationTheory& nt, int iter);

    /// Return the geometric mean volume [m3]
    double geometric_mean_v();

    /// Return the geometric standard deviation of the lognormal distribution
    double standard_dev() { return 3.0 * sqrt(ln2_standard_dev()); }

	/// Return the condenstion source term
	double get_cond_term() const { return M1_cond / species.m_volume(); }

	/// Get M1
	double get_M1()const { return M1; }

	/// Get M2
	double get_M2()const { return M2; }

	/// Get Lognormal
	void get_lognormal_val();

	/// Get Lognormal (Streamlines)
	void get_lognormal_val(int _s_index);

};


#endif // MOMENTMODELPRATSINIS_H
