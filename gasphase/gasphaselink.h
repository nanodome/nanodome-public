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

#ifndef GASPHASELINK_H
#define GASPHASELINK_H

#include "gasphase.h"
#include "../spline.h"

class GasPhaseLink : public GasPhase {

    /// Data structure describing the streamline imported from the CFD
    /// Spline describing the Temperature
    Spline *sT;

    /// Spline describing the pressure
    Spline *sP;

    /// Spline descibing the species molar fraction
    std::vector<Spline*> sC;

    /// Actual time for the gasphase
    double actual_time;

	/// Molar fraction
	double MF;

public:

    /// Provide the streamline with data for the GP evolution in time
    /// \param Spline* _T Spline describing the Temperature
    /// \param Spline* _P Spline describing the pressure
    /// \param Spline* _C Spline describing the species Concentration
    /// double s_time
    GasPhaseLink(Spline* _T, Spline* _P, std::vector<Species> _species,
                 std::vector<Spline*> _C, double s_time);

	/// Provide the streamline with data for the GP evolution in time
	/// \param Spline* _T Spline describing the Temperature
	/// \param Spline* _P Spline describing the pressure
	/// \param double s_time starting time of the streamline
	///	\param double nucl_species_start_m_frac starting value for the molar fraction provided by the streamline
	GasPhaseLink(Spline* _T, Spline* _P, std::vector<Species> _species, double s_time, 
				 double nucl_species_start_m_frac);

    /// Follow the streamline evolution using a timestep dt
    /// \param t simulation time [s]
	/// \param dt simulation timestep [s]
    void timestep(double t, double dt);

	/// Follow the streamline evolution using a timestep dt
	/// \param t simulation time [s]
	/// \param dt simulation timestep [s]
	void timestep_temp_grad(double t, double dt, std::valarray<double> w = { 0.0,0.0 });

	/// Get Molar fraction
	double get_molar_f() const { return MF; }



};

#endif // GASPHASELINK_H
