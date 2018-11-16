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

#include "momentmodelpratsinis.h"

#include <math.h>
#include <iostream>
#include <fstream>


MomentModelPratsinis::MomentModelPratsinis(Species _sp) : MomentModel(_sp) {

    M0=0; M1=0; M2=0;
}

MomentModelPratsinis::MomentModelPratsinis(Species _sp, double _M0, double _M1, double _M2) : MomentModel(_sp) {

	M0 = _M0; M1 = _M1; M2 = _M2;
}


double MomentModelPratsinis::timestep(double dt, const GasPhase& gp, const NucleationTheory& nt) {

    double T  = gp.get_T();
    double S  = gp.get_S("Si");

	/*if (S < 1.0)
		S = 1.0;*/

	double ns_test = gp.get_n("Si");

    double J = nt.nucleation_rate(T,S);
    double j = nt.stable_cluster_size(T,S);

    double gamma = gp.get_gamma();

    // temporary moments values
    double M0_, M1_, M2_;

    // volume of the critical size cluster
    // if S<1 the NucleationTheory object should return 1 as stable critical cluster size
    double vm = j * species.m_volume();

	/// dissipation limit
	double vm1 = 0.0;

    // nucleation source for each moment
    double M0_nucl = J;
    double M1_nucl = J * vm;
    double M2_nucl = J * vm*vm;

    double vg = geometric_mean_v();
    double ln2sg = ln2_standard_dev();

    // condensation source for each moment
	//double M1_cond;
	double M2_cond;
	M1_cond = csi1(T)*(S - 1)*M_k(2. / 3., ln2sg, vg);
	M2_cond = 2 * csi1(T)*(S - 1)*M_k(5. / 3., ln2sg, vg);

	// TEST
	/*M1_cond = 0.0;
	M2_cond = 0.0;*/
	

    double sg = exp(sqrt(ln2sg));

    // coagulation source for each moment
    double M0_coag =   csi2(T)*zeta0(sg)*(  M_k(2./3.,ln2sg,vg) * M_k(-1./2.,ln2sg,vg) +
                                          2*M_k(1./3.,ln2sg,vg) * M_k(-1./6.,ln2sg,vg) +
                                            M_k(1./6.,ln2sg,vg) * M0);

    double M2_coag = 2*csi2(T)*zeta2(sg)*(  M_k(5./3.,ln2sg,vg) * M_k(1./2.,ln2sg,vg) +
                                          2*M_k(4./3.,ln2sg,vg) * M_k(5./6.,ln2sg,vg) +
                                            M_k(7./6.,ln2sg,vg) * M1);

	// TEST
	M0_coag = 0.0;
	M2_coag = 0.0;

    // dissolution source: if S<1 then the flux of particles becoming smaller than a two monomers cluster
    // due to evaporation is removed from the particle set
    // in this case vm = 1 must hold
    double sigma = 3.0 * sqrt(ln2sg);
    double M0_diss = 0;
    double M1_diss = 0;
    double M2_diss = 0;

    if((S<1)&&(vg>0)&&(sigma>0)) {
		
		//vm *= 2;
		//
		//double nm = M0/(sigma*sqrt(2*M_PI)) * exp(-pow(log(vm/vg),2)/(2*sigma*sigma)) / vm;
  //      double G0 = csi1(T)*pow(vm,2./3.)*(S-1);

  //      M0_diss = -G0*nm; // particle number reduction
  //      M1_diss = -G0*nm*vm; // particle volume reduction
  //      M2_diss = -G0*nm*vm*vm; // M2 reduction

		vm1 = vg*std::exp(-1.82*1.41*sigma);

		double nm = M0 / (sigma*sqrt(2 * M_PI)) * exp(-pow(log(vm1 / vg), 2) / (2 * sigma*sigma)) / vm1;
		double G0 = csi1(T)*pow(vm1, 2. / 3.)*(S - 1);

		M0_diss = -G0*nm; // particle number reduction
		M1_diss = -G0*nm*vm1; // particle volume reduction
		M2_diss = -G0*nm*vm1*vm1; // M2 reduction



		////M0_diss *= 0;
		//M1_diss *= 1.0e5;
		//M2_diss *= 1.0e5;

    }

    // moments equations solution (simple first order explicit method) (TO BE IMPROVED!!!)
    M0_ = M0 + (M0_nucl           - M0_coag - M0_diss) * dt - gamma*M0*dt;
    M1_ = M1 + (M1_nucl + M1_cond           - M1_diss) * dt - gamma*M1*dt;
    M2_ = M2 + (M2_nucl + M2_cond + M2_coag - M2_diss) * dt - gamma*M2*dt;

	/*M0_ = M0 + (M0_nucl - M0_coag - M0_diss) * dt;
	M1_ = M1 + (M1_nucl + M1_cond - M1_diss) * dt;
	M2_ = M2 + (M2_nucl + M2_cond + M2_coag - M2_diss) * dt;*/

    // nucleating species consumption for
    // + homogeneous nucleation: J*j
    // + heterogenous nucleation: M1_cond/vm
    // - dissolution of evaporating particles reducing their size below vm
    double g = J*j + M1_cond/species.m_volume() - M1_diss/species.m_volume();

	// test
	g = 0.0;

	// DEBUG
	if (std::isnan(M0_) || std::isinf(M0_)) {
		std::cout
			<< " M0: " << M0 << std::endl
			<< " M1: " << M1 << std::endl
			<< " M2: " << M2 << std::endl
			<< " M0_nucl: " << M0_nucl << std::endl
			<< " M0_coag: " << M0_coag << std::endl
			<< " M0_diss: " << M0_diss << std::endl
			<< " csi2(T): " << csi2(T) << std::endl
			<< " zeta0(sg): " << zeta0(sg) << std::endl
			<< " sg: " << sg << std::endl
			<< " vg: " << vg << std::endl
			<< " ln2sg: " << ln2sg << std::endl
			<< " M_k(2. / 3., ln2sg, vg): " << M_k(2. / 3., ln2sg, vg) << std::endl
			<< " M_k(-1. / 2., ln2sg, vg): " << M_k(-1. / 2., ln2sg, vg) << std::endl
			<< " M_k(1. / 3., ln2sg, vg): " << M_k(1. / 3., ln2sg, vg) << std::endl
			<< " M_k(-1. / 6., ln2sg, vg): " << M_k(-1. / 6., ln2sg, vg) << std::endl
			<< " M_k(1. / 6., ln2sg, vg): " << M_k(1. / 6., ln2sg, vg) << std::endl
			<< std::endl;
		std::cout
			<< "J: " << nt.nucleation_rate(gp.get_T(), gp.get_S("Si")) << std::endl;
		system("PAUSE");
	}

   // TEST
	if (M0_<0.0 || M1_ < 0.0 || M2_ < 0.0){
		M1_ = 0.0;
		M2_ = 0.0;
		M0_ = 0.0;
	}

	// update moments
    M0 = M0_; M1 = M1_; M2 = M2_;
	

    //std::cout << M0 << ' ' << M1 << ' ' << M2 << ' ' << vm << ' '<< vg << ' ' << sigma << ' ' << xs << ' ' << M0_diss << ' ' << M1_diss << ' ' << M2_diss << std::endl;

    return g;
}

double MomentModelPratsinis::timestep(double dt, const GasPhase& gp, const NucleationTheory& nt, int iter) {

	double T = gp.get_T();
	double S = gp.get_S("Si");

	double ns_test = gp.get_n("Si");

	double J = nt.nucleation_rate(T, S);
	double j = nt.stable_cluster_size(T, S);

	double gamma = gp.get_gamma();

	// temporary moments values
	double M0_, M1_, M2_;

	// volume of the critical size cluster
	// if S<1 the NucleationTheory object should return 1 as stable critical cluster size
	double vm = j * species.m_volume();

	// nucleation source for each moment
	double M0_nucl = J;
	double M1_nucl = J * vm;
	double M2_nucl = J * vm*vm;

	double vg = geometric_mean_v();
	double ln2sg = ln2_standard_dev();

	// condensation source for each moment
	double M1_cond;
	double M2_cond;
	M1_cond = csi1(T)*(S - 1)*M_k(2. / 3., ln2sg, vg);
	M2_cond = 2 * csi1(T)*(S - 1)*M_k(5. / 3., ln2sg, vg);


	double sg = exp(sqrt(ln2sg));

	// coagulation source for each moment
	double M0_coag = csi2(T)*zeta0(sg)*(M_k(2. / 3., ln2sg, vg) * M_k(-1. / 2., ln2sg, vg) +
										2 * M_k(1. / 3., ln2sg, vg) * M_k(-1. / 6., ln2sg, vg) +
										M_k(1. / 6., ln2sg, vg) * M0);

	double M2_coag = 2 * csi2(T)*zeta2(sg)*(M_k(5. / 3., ln2sg, vg) * M_k(1. / 2., ln2sg, vg) +
											2 * M_k(4. / 3., ln2sg, vg) * M_k(5. / 6., ln2sg, vg) +
											M_k(7. / 6., ln2sg, vg) * M1);

	//if (M2_coag < 0.0 ) M2_coag = 0.0;

	// dissolution source: if S<1 then the flux of particles becoming smaller than a two monomers cluster
	// due to evaporation is removed from the particle set
	// in this case vm = 1 must hold
	double sigma = 3.0 * sqrt(ln2sg);
	double M0_diss = 0;
	double M1_diss = 0;
	double M2_diss = 0;

	if ((S<1) && (vg>0) && (sigma>0)) {

		vm *= 2;

		double nm = M0 / (sigma*sqrt(2 * M_PI)) * exp(-pow(log(vm / vg), 2) / (2 * sigma*sigma)) / vm;
		double G0 = csi1(T)*pow(vm, 2. / 3.)*(S - 1);

		M0_diss = -G0*nm; // particle number reduction
		M1_diss = -G0*nm*vm; // particle volume reduction
		M2_diss = -G0*nm*vm*vm; // M2 reduction
	}

	// moments equations solution (simple first order explicit method) (TO BE IMPROVED!!!)
	M0_ = M0 + (M0_nucl - M0_coag - M0_diss) * dt - gamma*M0*dt;
	M1_ = M1 + (M1_nucl + M1_cond - M1_diss) * dt - gamma*M1*dt;
	M2_ = M2 + (M2_nucl + M2_cond + M2_coag - M2_diss) * dt - gamma*M2*dt;

	// nucleating species consumption for
	// + homogeneous nucleation: J*j
	// + heterogenous nucleation: M1_cond/vm
	// - dissolution of evaporating particles reducing their size below vm
	double g = J*j + M1_cond / species.m_volume() - M1_diss / species.m_volume();

	// DEBUG

	/*if (std::isnan(M0_) || std::isinf(M0_) || iter > 18272595) {
		std::cout
			<< " iter "<<iter<< std::endl
			<< " M0: " << M0 << std::endl
			<< " M1: " << M1 << std::endl
			<< " M2: " << M2 << std::endl
			<< " M0_nucl: " << M0_nucl << std::endl
			<< " M0_coag: " << M0_coag << std::endl
			<< " M0_diss: " << M0_diss << std::endl
			<< " csi2(T): " << csi2(T) << std::endl
			<< " zeta0(sg): " << zeta0(sg) << std::endl
			<< " sg: " << sg << std::endl
			<< " vg: " << vg << std::endl
			<< " ln2sg: " << ln2sg << std::endl
			<< " M_k(2. / 3., ln2sg, vg): " << M_k(2. / 3., ln2sg, vg) << std::endl
			<< " M_k(-1. / 2., ln2sg, vg): " << M_k(-1. / 2., ln2sg, vg) << std::endl
			<< " M_k(1. / 3., ln2sg, vg): " << M_k(1. / 3., ln2sg, vg) << std::endl
			<< " M_k(-1. / 6., ln2sg, vg): " << M_k(-1. / 6., ln2sg, vg) << std::endl
			<< " M_k(1. / 6., ln2sg, vg): " << M_k(1. / 6., ln2sg, vg) << std::endl
			<< std::endl;
		std::cout
			<< " J: " << nt.nucleation_rate(gp.get_T(), gp.get_S("Si")) << std::endl;

		system("PAUSE");
	}*/

	// update moments
	M0 = M0_; M1 = M1_; M2 = M2_;


	//std::cout << M0 << ' ' << M1 << ' ' << M2 << ' ' << vm << ' '<< vg << ' ' << sigma << ' ' << xs << ' ' << M0_diss << ' ' << M1_diss << ' ' << M2_diss << std::endl;

	return g;
}


double MomentModelPratsinis::M_k(double k, double ln2sg, double vg) {

    return (vg>0.0) ? M0*pow(vg,k)*exp(4.5*k*k*ln2sg) : 0.0;
}


double MomentModelPratsinis::geometric_mean_v() {

    return ((M0>0.0)&&(M2>0.0)) ? M1*M1/(pow(M0,1.5)*sqrt(M2)) : 0.0;
}


double MomentModelPratsinis::ln2_standard_dev() {

    double value = 0.0;

    double arg = M0*M2/(M1*M1);

    if(arg>1.0)
        value = (1./9.) * log(arg);
	// Debug
	if (std::isinf(value)) {
		std::cout << "SG COMPUTATION" << std::endl;
		std::cout
			<< " M0: " << M0 << std::endl
			<< " M1: " << M1 << std::endl
			<< " M2: " << M2 << std::endl;
	}

    return value;
}


double MomentModelPratsinis::zeta0(double sg) {

    return 0.633 + 0.092*sg*sg - 0.022*sg*sg*sg;
}


double MomentModelPratsinis::zeta2(double sg) {

    return 0.39 + 0.5*sg - 0.214*sg*sg + 0.029*sg*sg*sg;
}


double MomentModelPratsinis::csi1(double T) {

    return species.m_volume()*species.n_sat(T)*4.835975862049408*sqrt(K_BOL*T/(2*M_PI*species.get_mass()));
}


double MomentModelPratsinis::csi2(double T) {

    return 0.787623317899743*sqrt(6*K_BOL*T/species.get_bulk_density(T));
}

void MomentModelPratsinis::get_lognormal_val() {

	std::ofstream _o_file;
	_o_file.open("Lognormal.dat", std::ofstream::app);

	// Moment 0
	(std::isnan(M0) || std::isinf(M0)) ? _o_file << 0.0 : _o_file << M0;
	_o_file << " ";

	// Moment 1
	(std::isnan(M1) || std::isinf(M1)) ? _o_file << 0.0 : _o_file << M1;
	_o_file << " ";

	// Moment 2
	(std::isnan(M2) || std::isinf(M2)) ? _o_file << 0.0 : _o_file << M2;
	_o_file << " ";

	//vg
	double vg = geometric_mean_v();
	(std::isnan(vg) || std::isinf(vg)) ? _o_file << 0.0 : _o_file << vg;
	_o_file << " ";

	/// sigma g
	double ln2sg = ln2_standard_dev();
	(std::isnan(ln2sg) || std::isinf(ln2sg)) ? _o_file << 0.0 : _o_file << ln2sg;
	_o_file << " ";
	

	_o_file << std::endl;

	// close file
	_o_file.close();

}

void MomentModelPratsinis::get_lognormal_val(int _s_index) {

	std::ofstream _o_file;
	_o_file.open("Lognormal_"+std::to_string(_s_index) + ".dat", std::ofstream::app);

	// Moment 0
	(std::isnan(M0) || std::isinf(M0)) ? _o_file << 0.0 : _o_file << M0;
	_o_file << " ";

	// Moment 1
	(std::isnan(M1) || std::isinf(M1)) ? _o_file << 0.0 : _o_file << M1;
	_o_file << " ";

	// Moment 2
	(std::isnan(M2) || std::isinf(M2)) ? _o_file << 0.0 : _o_file << M2;
	_o_file << " ";

	//vg
	double vg = geometric_mean_v();
	(std::isnan(vg) || std::isinf(vg)) ? _o_file << 0.0 : _o_file << vg;
	_o_file << " ";

	/// sigma g
	double ln2sg = ln2_standard_dev();
	(std::isnan(ln2sg) || std::isinf(ln2sg)) ? _o_file << 0.0 : _o_file << ln2sg;
	_o_file << " ";


	_o_file << std::endl;

	// close file
	_o_file.close();

}
