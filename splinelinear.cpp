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

#include "splinelinear.h"

SplineLinear::SplineLinear() :Spline() {

	alglib::spline1dbuildlinear(alg_x, alg_y, s);
}

SplineLinear::SplineLinear(std::vector<double> _x, std::vector<double> _y):
	Spline(_x, _y) {

	alglib::spline1dbuildlinear(alg_x, alg_y, s);

}

double SplineLinear::get_v(double _x) {

	if((x_min <= _x) && (_x <= x_max))
		return alglib::spline1dcalc(s, _x);
	else{
		i_o err;
		int std_prec = std::cout.precision();
		std::cout.precision(15);
		std::string error_entry = "Linear Spline requested value: " + std::to_string(_x) + "out of limits[" + std::to_string(x_min) + " " + std::to_string(x_max) + "]\n";
		std::cout.precision(std_prec);
		err.error_entry_not_blocking(error_entry);
	}

}
