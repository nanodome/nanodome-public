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

#include "momentmodelfriedlander.h"


MomentModelFriedlander::MomentModelFriedlander(Species _sp) : MomentModel(_sp) {

    M0=0.; M1=0.; M2=0.;
}


double MomentModelFriedlander::timestep(double dt,  const GasPhase& gp, const NucleationTheory& nt) {

    double T  = gp.get_T();
    double S  = gp.get_S("Si");

    double J = nt.nucleation_rate(T,S);
    double j = nt.stable_cluster_size(T,S);

    // volume exansion coefficient
    double gamma = gp.get_gamma();

    // temporary moments values
    double M0_, M1_, M2_;

    double s = species.m_surface(); //monomer area
    double v = species.m_volume(); // monomer volume
    double d = 2.*pow( (3./4.) * v * j / M_PI, 1./3.); // stable cluster diameter

    double b1 = 2. * species.n_sat(T) * v * sqrt(K_BOL*T/(2.*M_PI*species.get_mass()));

    // nucleation source for each moment
    double M0_nucl = J;
    double M1_nucl = J * d;
    double M2_nucl = J * pow(j,(2./3.)) * s;

    // condensation source for each moment
    double M1_cond = (S-1)*b1*M0;
    double M2_cond = 2*M_PI*b1*(S-1)*M1;

	// test
	M1_cond_test = M1_cond;

    // moments method equations (simple first order explicit method) (TO BE IMPROVED!!!)
    M0_ = M0 + (M0_nucl          ) * dt  - gamma*M0*dt;
    M1_ = M1 + (M1_nucl + M1_cond) * dt  - gamma*M1*dt;
    M2_ = M2 + (M2_nucl + M2_cond) * dt  - gamma*M2*dt;

    // nucleating species consumption for
    // - homogeneous nucleation: J*j
    // - heterogenous nucleation: (S-1.)*b1*A/(2.*v)
    double g = J*j + (S-1.)*b1*M2/(2.*v);

    // update moments
    M0 = M0_; M1 = M1_; M2 = M2_;

    return g;
}

