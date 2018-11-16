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

#ifndef FREEFLOWKERNEL_H
#define FREEFLOWKERNEL_H

#include "binarykernel.h"


template<typename T>
class FMBinaryKernel : public BinaryKernel<T> {

public:

    double get_coag_rate(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp);

    std::pair<std::shared_ptr<T>,std::shared_ptr<T>> get_coag_particles(std::list<std::shared_ptr<T>>& objects, const GasPhase& gp);
};

#include "fmbinarykernel.cpp"

#endif // FREEFLOWKERNEL_H
