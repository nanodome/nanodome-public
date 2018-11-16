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

#ifndef _SECTIONAL_H
#define _SECTIONAL_H

#include "../gasphase/gasphasecv.h"
#include "../cnt.h"
#include <valarray>

class SectionalModel {

	/// Number of bins of the method
	int bins;

	/// Scale factor
	double scale_f;

	/// Minimum feaseble volume
	double v_min;

	/// Volume of the bins [m3]
	//std::valarray<double> bin_values;
	double* bin_values;

	/// Diameters of the bins [m]
	//std::valarray<double> bin_diam;
	double* bin_diam;

	/// Cardinality of the bins
	//std::valarray<double> bin_card;
	double* bin_card;

	/// Beta factors look-up table
	double *beta_factors;

	/// chi factors lookup table
	double *chi_factors;

	/// Species(Monospecies) in the simulation
	Species species;

	/// TEST PURPOSES
	double gi_nucl;
	double gi_cond;
	double gi_cond_prakash;


	/// Nucleation
	double nucleation(double _dt, GasPhaseCV& _gp, ClassicalNucleationTheory& _cnt);

	/// Coagulation
	void coagulation(double _dt);

	/// Condensation
	double condensation(double _dt, GasPhaseCV& _gp);

	/// Csi parameter
	double chi(int _i, int _j, int _k);

	/// beta parameter (Friedlander Free molecular)(coagulation)
	double beta(int _i, int _j, int _T, double _si_density);

	/// beta parameter (Fuchs Free molecular + Continuum) (coagulation)
	double beta2(int _i, int _j, int _T, double _si_density, double _mfp, double _si_mass, double _mon_vol);

	/// beta parameter (Friedlander Brownian Motion) (coagulation)
	double beta3(int _i, int _j);

	/// beta condensation
	double beta_cond(int _k, double _T, double _mon_vol, double _si_density);


	/// csi parameter
	double csi(int _i, int _k, double _delta_vol, double& _delta_ik );

public:

	/// Constructor
	///	\param: int _bins: number of bins in the simulation
	///	\param: double _scale_f: bins scale factor
	///	\param double _v_min: minimum volume feaseble
	SectionalModel(int _bins, double _scale_f, double _v_min, Species _sp);

	/// Timestep
	double timestep(double dt, GasPhase& _gp, ClassicalNucleationTheory& _cnt, std::string _nuc_species);

	/// Density [#/m3]
	double get_density();

	/// Mean Diameter
	double get_mean_diameter();

	/// Kelvin Effect
	double Kelvin_effect(double _T, int _bin_idx);

	/// Condensation Prakash
	double cond_prakash(int _k, double _T, double _mon_vol, double _ns, double _S);

	/// Gasphase consuption Prakash
	double cons_cond_prakash(int _k, double _T, double _mon_vol, double _ns);

	/// print bin cardinalities
	void print_bins(double _t);

	/// print bins cardinalities (Streamlines)
	void print_bins(double _t, int _stream_idx);

	/// print bin values
	void print_bins_values();

	/// print bin values (Streamline)
	void print_bins_values(int _stream_idx);

	/// print consuptions
	double get_gi_cond()const { return gi_cond; 
		//return gi_cond_prakash; 
	}
	double get_gi_nucl()const { return gi_nucl; }

	/// Create dummy distribution
	void dummy_distribution(double _vg, double _ln2sg, double _M0);
};





#endif
