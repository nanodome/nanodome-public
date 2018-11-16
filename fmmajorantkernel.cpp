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

#include "fmmajorantkernel.h"
#include "ndm_random.h"


template<typename T>
double FMMajorantKernel<T>::get_coag_rate(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp) {

    double d_col;
    double mass;

    K1 = K2 = K21 = K22 = 0;

    for(auto& a: objects) {

        d_col = a->get_collision_diameter();
        mass = a->get_mass();

        d_col = square(d_col);
        mass = 1./sqrt(mass);

        K1 += d_col*mass;
        K21 += d_col;
        K22 += mass;
    }

    K1 *= 2.0 * (objects.size() - 2.0);
    K2 = 2.0 * K21 * K22;

    double k_maj = 2;

    return k_maj * sqrt(0.5*M_PI*K_BOL*gp.get_T())*(K1+K2);
}


template<typename T>
std::pair<std::shared_ptr<T>,std::shared_ptr<T>> FMMajorantKernel<T>::get_coag_particles(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp) {

    std::pair<std::shared_ptr<T>,std::shared_ptr<T>> ij;

    double rho = ndm::uniform_double_distr(ndm::rand_gen);

    if(rho<K1/(K1+K2)) {

        while(ij.first==ij.second) {

            // pick a random a0 with w = 1
            auto it = objects.begin();

            double w0_rnd = ndm::uniform_double_distr(ndm::rand_gen);

            ij.first = *(std::next(it,int(w0_rnd*(objects.size()-1))));

            // pick a random a1 with w = d^2/sqrt(m) using Cumulative Density Func
            double w1 = 0;
            double w1_rnd = ndm::uniform_double_distr(ndm::rand_gen);
            double w1_norm = 0.5*K1/(objects.size()-2.0);

            for(auto it = objects.begin(); it!=objects.end(); ++it) {

                double d_col = (*it)->get_collision_diameter();
                double mass = (*it)->get_mass();

                d_col = square(d_col); mass = 1./sqrt(mass);
                w1 += d_col*mass/w1_norm;

                if(w1>w1_rnd) {
                    ij.second = *(it);
                    break;
                }
            }
        }

    } else {

        while(ij.first==ij.second) {

            double w0 = 0;
            double w1 = 0;

            double w0_rnd = ndm::uniform_double_distr(ndm::rand_gen);
            double w1_rnd = ndm::uniform_double_distr(ndm::rand_gen);

            double w0_norm = K21;
            double w1_norm = K22;

            // pick a random a0 with w = d^2
            for(auto it=objects.begin(); it!=objects.end(); ++it) {

                double d_col = (*it)->get_collision_diameter();

                d_col = square(d_col);
                w0 += d_col/w0_norm;

                if(w0>w0_rnd) {
                    ij.first = *(it);
                    break;
                }
            }

            // pick a random a1 with w = 1/sqrt(m)
            for(auto it=objects.begin(); it!=objects.end(); ++it) {

                double mass = (*it)->get_mass();

                mass = 1./sqrt(mass);
                w1 += mass/w1_norm;

                if(w1>w1_rnd) {
                    ij.second = *(it);
                    break;
                }
            }
        }
    }

    // calculate kernels
    double m0 = ij.first->get_mass();
    double m1 = ij.second->get_mass();
    double d0 = ij.first->get_collision_diameter();
    double d1 = ij.second->get_collision_diameter();

    // free molecular kernel
    double K = sqrt(1/m0 + 1/m1) * square(d0 + d1);

    // free molecular majorant kernel
    double k_maj = 2.0;
    double K_maj = k_maj * (1/sqrt(m0) + 1/sqrt(m1)) * (square(d0) + square(d1));

    // check if it is a real or a fictious coagulation
    double delta = ndm::uniform_double_distr(ndm::rand_gen);
    if(delta > K/K_maj) {
        ij.first.reset();
        ij.second.reset();
    }

    return ij;
}
