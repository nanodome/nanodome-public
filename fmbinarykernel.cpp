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

#include "fmbinarykernel.h"
#include "collisionkernel.h"
#include "ndm_random.h"


template<typename T>
double FMBinaryKernel<T>::get_coag_rate(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp) {

    double R = 0.0;

    // sum all the rates for each possible binary interaction
    for(auto i = objects.begin(); i!=objects.end(); ++i)
        for(auto j = std::next(i); j!=objects.end(); ++j)
            R += FreeMolecularKernel(**i,**j,gp.get_T());

    return R;
}


template<typename T>
std::pair<std::shared_ptr<T>,std::shared_ptr<T>> get_coag_particles(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp) {

    std::pair<std::shared_ptr<T>,std::shared_ptr<T>> ij;

    // recalculate total coagulation rate
    double R = get_coag_rate(objects,gp);

    double rho = ndm::uniform_double_distr(ndm::rand_gen);

    // pick a random couple i,j with w = Rij using Cumulative Density Funct
    double R_cum = 0;

    for(auto i = objects.begin(); i!=objects.end(); ++i) {
        for(auto j = std::next(i); j!=objects.end(); ++j) {

            double Rij = FreeMolecularKernel(**i,**j,gp.get_T());

            R_cum += Rij/R;

            if(R_cum>rho) {
                ij.first = *i;
                ij.second = *j;
                return ij;
            }
        }
    }
}
