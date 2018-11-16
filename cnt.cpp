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

#include "cnt.h"

#include <cmath>
#include <limits>
#include <exception>
#include <stdexcept>
#include <iostream>

ClassicalNucleationTheory::ClassicalNucleationTheory(Species _species) :
    species(_species) { }


double ClassicalNucleationTheory::nucleation_rate(double T, double S)  const {

    double rate = 0.0;
	//unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

    // check if species is saturated, if not nucleation rate is left to zero
    if(S>1.0) {

        // normalized surface tension
        double theta = species.s_ten(T)*species.m_surface()/(K_BOL*T);
        double ns_sat = species.n_sat(T);

        double A = (S*ns_sat*ns_sat);
        double B = species.m_volume();
        double C1 = 2.0*species.s_ten(T);
        double C2 = M_PI*species.get_mass();
        double D1 = theta - 4.0*pow(theta, 3) / (27.0*pow(log(S), 2));
        double C = sqrt(C1 / C2);

		double D = exp(D1);
		
		rate = A * B * C * D;

    } else if (S==1) {

        rate = std::numeric_limits<double>::min();
    }

    return rate;
}


double ClassicalNucleationTheory::stable_cluster_size(double T, double S) const {

    double c_size = 1.0;

    // check if species is saturated; if not, the stable size has no sense and is set to one
    // meaning that the smallest cluster is a single monomer (no-cluster)
    if(S>1) {
        c_size = 2.0 * species.m_surface() * species.s_ten(T) / (3*K_BOL*T*log(S));
        c_size = pow(c_size,3);
		}

    return c_size;
}


double ClassicalNucleationTheory::stable_cluster_diameter(double T, double S) const {

    return 2*pow( (3./4.) * species.m_volume() * stable_cluster_size(T,S) / M_PI, 1./3.);
}


double ClassicalNucleationTheory::condensation_rate(double T, double S) const {

    return species.p_sat(T)*(S-1.0) / sqrt(2*M_PI*species.get_mass()*K_BOL*T);
}
