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

#include "sectional.h"

#include <iostream>
#include <fstream>

SectionalModel::SectionalModel(int _bins, double _scale_f, double _v_min, Species _sp):
	bins(_bins), scale_f(_scale_f), v_min(_v_min), species(_sp){

	/*bin_values.resize(bins);
	bin_diam.resize(bins);*/
	bin_values = (double*)calloc(bins, sizeof(double));
	bin_diam = (double*)calloc(bins, sizeof(double));

	// Calculate bin values and diameters
	bin_values[0] = v_min;
	bin_diam[0] = std::pow((6.0*v_min) / M_PI, 1.0 / 3.0);

	for (size_t b = 1; b < bins; b++) {
		bin_values[b] = bin_values[b - 1] * scale_f;
		double bin_diam_v = std::pow((6.0*bin_values[b]) / M_PI, 1.0 / 3.0);
		bin_diam[b] = bin_diam_v;
	}

	//bin_card.resize(bins, 0.0);
	bin_card = (double*)calloc(bins, sizeof(double));

	/// Lookup table creation for beta (coagulation)
	beta_factors = (double*)malloc(bins*bins * sizeof(double));
	for (int i = 0; i < bins; i++) {
		for (int j = 0; j < bins; j++) {
			double F3 = sqrt(1.0 / bin_values[i] + 1.0 / bin_values[j]);
			double F4 = square(std::pow(bin_values[i], 1.0 / 3.0) + std::pow(bin_values[j], 1.0 / 3.0));
			beta_factors[i*bins + j] = F3*F4;
		}
	}

	/// Lookup table for chi factor (coagulation)
	chi_factors = (double*)malloc(bins*bins*bins * sizeof(double));

	for (int i = 0; i < bins; i++) {
		for (int j = 0; j < bins; j++) {
			for (int k = 0; k < bins; k++) {
				chi_factors[i*bins*bins + j*bins + k] = chi(i, j, k);
			}
		}
	}

	// TEST!!!
	//bin_card[30] = 1.0e18;
	//bin_card[19] = 1000;

}

double SectionalModel::timestep(double _dt, GasPhase& _gp, ClassicalNucleationTheory& _cnt, std::string _nuc_species) {

	double gi = 0.0; //[# monomers]

	//std::valarray<double> bin_source_nucl(0.0, bins); //[# particles]
	//std::valarray<double> bin_source_cond(0.0, bins); //[# particles]
	//std::valarray<double> bin_source_coag(0.0, bins); //[# particles]

	double* bin_source_nucl = (double*)calloc(bins, sizeof(double));
	double* bin_source_cond = (double*)calloc(bins, sizeof(double));
	double* bin_source_coag = (double*)calloc(bins, sizeof(double));
													  
	// nucleation
	double T = _gp.get_T();
	double S = _gp.get_S(_nuc_species);
	double ns = _gp.get_n(_nuc_species);
	
	/*if (S < 1.0) 
		S = 1.0;*/

	double J = _cnt.nucleation_rate(T, S);
	double j = _cnt.stable_cluster_size(T, S);

	// Elimination of particles smaller than the actual stable cluster size
	double vm = j * species.m_volume(); //[m3]

	/// Identify bin (USE BINARY SEARCH!!!)
	
	for (size_t b = 0; b < bins; b++) {
		// bin selected
		if (vm <= bin_values[b]) {
			// calculate weight
			double csi = 0.0;
			csi = vm / bin_values[b]; 
			// Update bin
			bin_source_nucl[b] += csi*J*_dt;
			break;
		}
	}
	
	// condensation

	
	gi_cond_prakash = 0.0;
	
	/// V 3.0
	/*
	double monomer_vol = species.m_volume();

	for (size_t k = 0; k < bins; k++) {
		double delta_p = cond_prakash(k, T, monomer_vol, ns, S);
		bin_source_cond[k] += delta_p*_dt;
		//gi_cond += delta_p*_dt*bin_values[k];
		gi_cond_prakash -= cons_cond_prakash(k, T, monomer_vol, ns);
	}
	
	*/
	
	// V 2.0
	gi_cond = 0.0;

	double Fs = _cnt.condensation_rate(T, S); //[# mon/m2 sec]

	double mon_flux = Fs*species.m_volume(); // [m3/m2 sec]
	
	for (size_t k = 0; k < bins; k++) {
		for (size_t i = 0; i < bins; i++) {
			double bin_v = bin_values[i]; //[m3] bin volume
			double bin_surf = 4.0*M_PI*pow((3.0*bin_v) / (4.0*M_PI), 2.0 / 3.0); // [m2] bin surface
			double delta_vol = mon_flux*bin_surf; // [m3/s]
			double bin_update = bin_v + delta_vol*_dt;
			gi_cond += delta_vol*_dt*bin_card[i]; // [m3] volume of condensating gas in dt
			double delta_ik;

			(k != i) ? delta_ik = 0.0 : delta_ik = 1.0;
			
			double csi_val = csi(i, k, delta_vol*_dt, delta_ik);

			//if (bin_update == bin_values[i] && i == k)
				//csi_val = 0.0; delta_ik = 0.0;

			
			bin_source_cond[k] += (bin_card[i] * (csi_val - delta_ik));
			
			//bin_source_cond[k] += (bin_card[i] * (csi(i, k, delta_vol) - delta_ik)) / _dt;
		}
	}
	

	/*
	/// V. 1.0
	for (int b = 0; b < bins; b++) {
		if (bin_card[b] > 1.0) {

			bin_source_cond[b] = bin_source_cond[b] - bin_card[b];

			double bin_v = bin_values[b]; //[m3] bin volume
			double bin_surf = 4.0*M_PI*pow((3.0*bin_v) / (4.0*M_PI), 2.0 / 3.0); // [m2] bin surface
			double delta_vol = mon_flux*bin_surf; // [m3/s]

			gi_cond += delta_vol*_dt*bin_card[b]; // [m3] volume of condensating gas in dt

			double bin_update = bin_v + delta_vol * _dt; //[m3] new volume for the bin with surface condensation in dt

			if (bin_update >= bin_values[bins - 1]) {
				bin_source_cond[bins - 1] = (bin_update*bin_card[b]) / bin_values[bins - 1]*_dt;
			}
			else {
				/// Search the new bin (Binary search)
				for (int n_b = 0; n_b < bins; n_b++) {
					if (bin_update < bin_values[n_b]) {
						// calculate weight
						double csi = (bin_values[n_b] - bin_update) / (bin_values[n_b] - bin_values[n_b - 1]);
						bin_source_cond[n_b - 1] += (1 - csi)*bin_card[b]*_dt;
						bin_source_cond[n_b] += csi*bin_card[b]*_dt;
						break;
					}

				}
			}
		}
	}
	*/
	// Nucleation Consumption

	// number of monomers nucleated from the gas phase
	for (size_t b = 0; b < bins; b++) {
		gi += (bin_source_nucl[b] * bin_values[b] / species.m_volume()) / _dt; // [# monomers / s]
	}

	gi_nucl=gi;

	// Condensation Consumption
	// number of monomers condensated from the gas phase due to surface condensation
	gi_cond = (gi_cond / species.m_volume()) / _dt;
	
	gi += gi_cond; //[# monomers / sec] 
	//gi += gi_cond_prakash;

	
	
	// coagulation
	double si_density = species.get_bulk_density(T);
	double mfp = _gp.get_mfp();
	double si_mass = species.get_mass();
	double si_vol = species.m_volume();

	/*
	// V 1.0
	for (size_t b1 = 0; b1 < bins; b1++) {
		for (size_t b2 = b1; b2 < bins; b2++) {

			if (bin_card[b1] != 0.0 && bin_card[b2] != 0.0) {

				double F1 = std::pow(3.0 / (4.0*M_PI), 1.0 / 6.0);
				double F2 = sqrt((6.0*K_BOL*T) / si_density);
				double F3 = sqrt(1.0 / bin_values[b1] + 1.0 / bin_values[b2]);
				double F4 = square(std::pow(bin_values[b1], 1.0 / 3.0) + std::pow(bin_values[b2], 1.0 / 3.0));
				double beta = F1*F2*F3*F4*_dt;
				// test

				double coag_vol = bin_values[b1] + bin_values[b2];

				if (coag_vol >= bin_values[bins - 1]) {
					bin_source_coag[bins - 1] = beta * bin_card[b1] * bin_card[b2] * coag_vol / bin_values[bins - 1];

					bin_source_coag[b1] -= beta * bin_card[b1] * bin_card[b2];
					bin_source_coag[b2] -= beta * bin_card[b1] * bin_card[b2];
				}
				else {
					// BINARY SEARCH!!!!
					for (size_t n_b = 0; n_b < bins; n_b++) {
						if (coag_vol < bin_values[n_b]) {

							double chi = (bin_values[n_b] - coag_vol) / (bin_values[n_b] - bin_values[n_b - 1]);

							if (chi > 1.0)
								std::cout << "MUORI ANCORA\n";

							bin_source_coag[n_b - 1] += (1 - chi) * beta * bin_card[b1] * bin_card[b2];
							bin_source_coag[n_b] += chi * beta * bin_card[b1] * bin_card[b2];

							if (std::isnan(bin_source_coag[n_b - 1]) || std::isinf(bin_source_coag[n_b - 1]) ||
								std::isnan(bin_source_coag[n_b]) || std::isinf(bin_source_coag[n_b])) {
								std::cout << "MUORI" << std::endl;
							}

							bin_source_coag[b1] -= beta * bin_card[b1] * bin_card[b2];
							bin_source_coag[b2] -= beta * bin_card[b1] * bin_card[b2];
							break;
						}
					}
				}
			}
		}
	}
	*/

	/*
	// V2.0
	for (size_t k = 0; k < bins; k++) {
		double Nk = bin_card[k];
		for (size_t i = 0; i < bins; i++) {
			double Ni = bin_card[i];
			//double beta_ik = beta(i, k, T, si_density);
			double beta_ik = beta2(i, k, T, si_density, mfp, si_mass, si_vol);
			//double beta_ik = beta3(i, k);
			bin_source_coag[k] -= beta_ik*Ni*Nk*_dt;
			if (k > 0) {
				for (size_t j = 0; j < bins; j++) {
					double Nj = bin_card[j];
					//double beta_ij = beta(i, j, T, si_density);
					double beta_ij = beta2(i, j, T, si_density, mfp, si_mass, si_vol);
					//double beta_ij = beta3(i, j);
					//bin_source_coag[k] += 0.5*beta_ij*chi(i, j, k)*Ni*Nj*_dt;
					bin_source_coag[k] += 0.5*beta_ij*chi_factors[i*bins*bins + j*bins + k]*Ni*Nj*_dt;
				}
			}
		}
	}
	*/

//#pragma omp parallel for collapse(2)
	//for (size_t k = 0; k < bins; k++) {
	//	for (size_t i = 0; i < bins; i++) {
	//		double Nk = bin_card[k];
	//		double Ni = bin_card[i];
	//		//double beta_ik = beta(i, k, T, si_density);
	//		double beta_ik = beta2(i, k, T, si_density, mfp, si_mass, si_vol);
	//		//double beta_ik = beta3(i, k);
	//		bin_source_coag[k] -= beta_ik*Ni*Nk*_dt;
	//		if (k > 0) {
	//			for (size_t j = 0; j < bins; j++) {
	//				double Nj = bin_card[j];
	//				//double beta_ij = beta(i, j, T, si_density);
	//				double beta_ij = beta2(i, j, T, si_density, mfp, si_mass, si_vol);
	//				//double beta_ij = beta3(i, j);
	//				//bin_source_coag[k] += 0.5*beta_ij*chi(i, j, k)*Ni*Nj*_dt;
	//				bin_source_coag[k] += 0.5*beta_ij*chi_factors[i*bins*bins + j*bins + k]*Ni*Nj*_dt;
	//			}
	//		}
	//	}
	//}
	
	
	// Update bins cardinality
	for (size_t b = 0; b < bins; b++) {
		
		double update_val = bin_card[b] + (bin_source_nucl[b] + bin_source_cond[b] + bin_source_coag[b]) 
						  - _gp.get_gamma()*bin_card[b]*_dt;

		//double update_val = bin_card[b] + (bin_source_nucl[b] + bin_source_cond[b] + bin_source_coag[b]);

		bin_card[b] = update_val;
		if (bin_card[b] < 0.0) bin_card[b] = 0.0;
	}

	// Test
	gi = 0.0;

	// TEST
	/*std::ofstream _o_file;
	_o_file.open("BinSource.dat", std::ofstream::app);

	for (size_t b = 0; b < bins; b++) {
		_o_file << bin_source_cond[b] << " ";
	}
	_o_file << std::endl;

	_o_file.close();*/

	return gi;
}

double SectionalModel::nucleation(double _dt, GasPhaseCV& _gp, 
								  ClassicalNucleationTheory& _cnt) {

	return 0.0;
}

double SectionalModel::condensation(double _dt, GasPhaseCV& _gp) {

	double gi = 0.0;
	

	return gi;

}

void SectionalModel::coagulation(double _dt) {



}

double SectionalModel::get_density() {

	double density = 0.0;
	for (size_t b = 0; b < bins; b++) {
		density += bin_card[b];
	}

	return density;
}

double SectionalModel::get_mean_diameter() {

	double mean_volume = 0.0;
	double tot_weights = 0.0;
	for (size_t b = 0; b < bins; b++) {
		if (bin_card[b] > 1.0) {
			tot_weights += bin_card[b];
			mean_volume += bin_card[b] * bin_values[b];
		}
	}

	if (tot_weights > 0) {
		mean_volume /= tot_weights;

		double mean_diameter = std::pow((3.0*mean_volume) / (4.0*M_PI), 1.0 / 3.0)*2.0;

		return mean_diameter;
	}
	else
		return 0.0;
}

void SectionalModel::print_bins(double _t) {

	std::ofstream _o_file;
	_o_file.open("BinValues.dat", std::ofstream::app);
	_o_file << _t << " ";
	for (size_t b = 0; b < bins; b++) {
		_o_file << bin_card[b] << " ";
	}
	_o_file << std::endl;

	_o_file.close();

}

void SectionalModel::print_bins(double _t, int _stream_idx) {

	std::ofstream _o_file;
	_o_file.open("BinValues_"+ std::to_string(_stream_idx) + "_Streamline.dat", std::ofstream::app);
	_o_file << _t << " ";
	for (size_t b = 0; b < bins; b++) {
		_o_file << bin_card[b] << " ";
	}
	_o_file << std::endl;

	_o_file.close();

}

void SectionalModel::print_bins_values() {

	std::ofstream _o_file;
	_o_file.open("BinValues.dat", std::ofstream::app);

	for (size_t b = 0; b < bins; b++) {
		_o_file << bin_values[b] << " ";
	}
	_o_file << std::endl;

	_o_file.close();

}

void SectionalModel::print_bins_values(int _stream_idx) {

	std::ofstream _o_file;
	_o_file.open("BinValues_" + std::to_string(_stream_idx) + "_Streamline.dat", std::ofstream::app);

	for (size_t b = 0; b < bins; b++) {
		_o_file << bin_values[b] << " ";
	}
	_o_file << std::endl;

	_o_file.close();

}

double SectionalModel::chi(int _i, int _j, int _k) {

	double coag_vol = bin_values[_i] + bin_values[_j];

	if (_k < bins - 1 || coag_vol < bin_values[_k]) {
		if (coag_vol >= bin_values[_k] && coag_vol <= bin_values[_k + 1]) {
			return (bin_values[_k + 1] - coag_vol) / (bin_values[_k + 1] - bin_values[_k]);
		}
		else if (coag_vol >= bin_values[_k - 1] && coag_vol <= bin_values[_k])
			return (coag_vol - bin_values[_k - 1]) / (bin_values[_k] - bin_values[_k - 1]);
		else
			return 0.0;
	}
	else {
		return 1.0;
	}
}

double SectionalModel::beta(int _i, int _j, int _T, double _si_density) {

	double F1 = 0.7876;//std::pow(3.0 / (4.0*M_PI), 1.0 / 6.0);
	double F2 = sqrt((6.0*K_BOL*_T) / _si_density);
	/*double F3 = sqrt(1.0 / bin_values[_i] + 1.0 / bin_values[_j]);
	double F4 = square(std::pow(bin_values[_i], 1.0 / 3.0) + std::pow(bin_values[_j], 1.0 / 3.0));*/
	double F3_F4 = beta_factors[_i*bins + _j];
	double beta = F1*F2*F3_F4;

	return beta;
}

double SectionalModel::beta_cond(int _k, double _T, double _mon_vol, double _si_density) {

	double F1 = std::pow(3.0 / (4.0*M_PI), 1.0 / 6.0);
	double F2 = sqrt((6.0*K_BOL*_T) / _si_density);
	double F3 = sqrt(1.0 / bin_values[_k] + 1.0 / _mon_vol);
	double F4 = square(std::pow(bin_values[_k], 1.0 / 3.0) + std::pow(_mon_vol, 1.0 / 3.0));
	double beta = F1*F2*F3*F4;

	return beta;

}

double SectionalModel::beta2(int _i, int _j, int _T, double _si_density, double _mfp, double _si_mass, double _mon_vol) {

	double beta_f = 0.0;

	// Particles diameters
	/*double di = 2.0*std::pow((3.0*bin_values[_i]) / (4.0*M_PI), 1.0 / 3.0);
	double dj = 2.0*std::pow((3.0*bin_values[_j]) / (4.0*M_PI), 1.0 / 3.0);*/
	double di = bin_diam[_i];
	double dj = bin_diam[_j];
	// Particles Knudsen Number
	double Ki = _mfp / (di*0.5);
	double Kj = _mfp / (dj*0.5);



	// Particles Diffusion coefficient
	double Fi1 = (K_BOL*_T) / (3.0*M_PI*ND_VISCOSITY*di);
	double Fi2 = 1 + Ki*(1.257 + 0.4*std::exp((-2.0*0.55) / Ki));
	double Dpi = Fi1 * Fi2;

	// Other factors
	double Fj1 = (K_BOL*_T) / (3.0*M_PI*ND_VISCOSITY*dj);
	double Fj2 = 1 + Kj*(1.257 + 0.4*std::exp((-2.0*0.55) / Kj));
	double Dpj = Fj1 * Fj2;

	/*

	double mi = _si_density*bin_values[_i];
	double mj = _si_density*bin_values[_j];

	*/
	
	double mon_vol = _mon_vol;

	double mi = bin_values[_i] / mon_vol*_si_mass;
	double mj = bin_values[_j] / mon_vol*_si_mass;
	

	double ci = std::sqrt((8.0*K_BOL*_T) / (M_PI*mi));
	double cj = std::sqrt((8.0*K_BOL*_T) / (M_PI*mj));

	double li = ((8.0*Dpi) / (M_PI*ci));
	double lj = ((8.0*Dpj) / (M_PI*cj));

	double gi = 1.0 / (3.0*di*li)*(std::pow(di + li, 3.0) - std::pow(square(di) + square(li), 3.0 / 2.0)) - di;
	double gj = 1.0 / (3.0*dj*lj)*(std::pow(dj + lj, 3.0) - std::pow(square(dj) + square(lj), 3.0 / 2.0)) - dj;

	// Beta calculation
	double G1 = 2.0*M_PI*(Dpi + Dpj)*(di + dj);
	double G2 = (di + dj) / ((di + dj) + 2.0*std::sqrt( square(gi) + square(gj) ));
	double G3 = (8.0*(Dpi + Dpj)) / ((di + dj)*std::sqrt( square(ci) + square(cj) ));

	beta_f = G1*1.0 / (G2 + G3);

	return beta_f;
}

double SectionalModel::beta3(int _i, int _j) {

	double beta = 0.0;

	beta = ((2.0*K_BOL) / (3.0*ND_VISCOSITY))*
		(1.0 / std::pow(bin_values[_i], 1.0 / 3.0) + 1.0 / std::pow(bin_values[_j], 1.0 / 3.0))*
		(std::pow(bin_values[_i], 1.0 / 3.0) + std::pow(bin_values[_j], 1.0 / 3.0));

	return beta;
}

double SectionalModel::csi(int _i, int _k, double _delta_vol, double& _delta_ik) {

	double csi = 0.0;

	if (_k < bins - 1 || bin_values[_i] + _delta_vol < bin_values[_k]) {

		if (bin_values[_k] <= bin_values[_i] + _delta_vol && bin_values[_k + 1] > bin_values[_i] + _delta_vol) {
			return csi = (bin_values[_k + 1] - (bin_values[_i] + _delta_vol)) / (bin_values[_k + 1] - bin_values[_k]);
		}
		else if (bin_values[_k - 1] < bin_values[_i] + _delta_vol && bin_values[_k] >= bin_values[_i] + _delta_vol) {
			return csi = (bin_values[_i] + _delta_vol - bin_values[_k - 1]) / (bin_values[_k] - bin_values[_k - 1]);
		}
		else
			return csi = 0.0;
	}
	else {
		if (_i != _k) {
			return csi = 1.0;
			//return csi = (bin_values[_i] + _delta_vol) / bin_values[_k];
		}
		else {
			_delta_ik = 0.0;
			return csi = 0.0; // volume
			//return csi = (bin_values[_i] + _delta_vol) / bin_values[_k];
		}
	}
	
}

double SectionalModel::Kelvin_effect(double _T, int _bin_idx) {

	double K_effect = 0.0;

	K_effect = species.n_sat(_T)*std::exp((4.0*species.s_ten(_T)*species.m_volume()) / (bin_diam[_bin_idx] * K_BOL*_T));

	return K_effect;

}

double SectionalModel::cond_prakash(int _k, double _T, double _mon_vol, double _ns, double _S) {

	double delta_part = 0.0;

	double Nk_1s = Kelvin_effect(_T, _k - 1);
	double Nks = Kelvin_effect(_T, _k);
	double Nkplus1 = Kelvin_effect(_T, _k + 1);

	double beta_k_1 = beta_cond(_k - 1, _T, _mon_vol, species.get_bulk_density(_T));
	double beta_k = beta_cond(_k, _T, _mon_vol, species.get_bulk_density(_T));
	double beta_kplus1 = beta_cond(_k + 1, _T, _mon_vol, species.get_bulk_density(_T));

	// Condensation

	if (_ns > Nk_1s && _k > 0)
		delta_part += (_mon_vol / (bin_values[_k] - bin_values[_k - 1]))*beta_k_1*(_ns - Nk_1s)*bin_card[_k - 1];
	if (_ns > Nks && _k < bins - 1)
		delta_part += -(_mon_vol / (bin_values[_k + 1] - bin_values[_k]))*beta_k*(_ns - Nks)*bin_card[_k];

	// Evaporation

	if (_ns < Nkplus1 && _k < bins - 1)
		delta_part += -(_mon_vol / (bin_values[_k + 1] - bin_values[_k]))*beta_kplus1*(_ns - Nkplus1)*bin_card[_k + 1];
	if (_ns < Nks && _k > 0)
		delta_part += (_mon_vol / (bin_values[_k] - bin_values[_k - 1]))*beta_k*(_ns - Nks)*bin_card[_k];
	if (_ns < Nks && _k == 0) {

		delta_part += (_mon_vol / (bin_values[_k] - _mon_vol))*beta_k*(_ns - Nks)*bin_card[_k];
	}

	return delta_part;


}

double SectionalModel::cons_cond_prakash(int _k, double _T, double _mon_vol, double _ns) {

	double gk = 0.0;

	double Nk_1s = Kelvin_effect(_T, _k - 1);
	double Nks = Kelvin_effect(_T, _k);
	double Nkplus1 = Kelvin_effect(_T, _k + 1);

	//test

	double beta_k_1 = beta_cond(_k - 1, _T, _mon_vol, species.get_bulk_density(_T));
	double beta_k = beta_cond(_k, _T, _mon_vol, species.get_bulk_density(_T));
	double beta_kplus1 = beta_cond(_k + 1, _T, _mon_vol, species.get_bulk_density(_T));

	// Condensation
	if (_ns > Nk_1s && _k > 0)
		gk += -beta_k_1*(_ns - Nk_1s)*bin_card[_k - 1];
	if (_ns > Nks)
		gk += -beta_k*(_ns - Nks)*bin_card[_k];

	// Evaporation
	if(_ns < Nkplus1 && _k < bins - 1)
		gk += -beta_kplus1*(_ns - Nkplus1)*bin_card[_k + 1];
	if (_ns < Nks)
		gk += -beta_k*(_ns - Nks)*bin_card[_k];


	return gk;

}

void SectionalModel::dummy_distribution(double _vg, double _ln2sg, double _M0) {

	for (size_t b = 1; b < bins; b++) {
		double F1 = 1.0 / (std::sqrt(2.0*M_PI*9.0*_ln2sg)*bin_values[b]);
		double num = square(std::log(bin_values[b] / _vg));
		double F2 = std::exp(-(num) / (18.0*_ln2sg));
		bin_card[b] = _M0*F1*F2*(bin_values[b] - bin_values[b-1]);
	}

}
