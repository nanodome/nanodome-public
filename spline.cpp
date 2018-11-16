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

#include "spline.h"

void Spline::scrubbing() {

	int n = x.size();

	// Scrubbing 
	for (int i = 1; i <= n-1; i++) {
		if ((x[i] - x[i - 1]) <= 0.0) {
			std::cout << "Scrubbing x: " << x[i] << "value" << std::endl;
			x.erase(x.begin() + (i-1));
			y.erase(y.begin() + (i-1));
			n = x.size();
		}
	}

	nodes = x.size();
}

double* convert_data(std::vector<double> _in) {

	double *res;
	int D = _in.size();

	res = (double*)malloc(D * sizeof(double));

	for (int i = 0; i < D; i++) {
		res[i] = _in[i];
	}

	return res;
}

Spline::Spline() {

	x_max = x_min = 0.0;
	alg_x = "[0, 1]" ;
	alg_y = "[0, 0]" ;
	nodes = 0;
}


Spline::Spline(std::vector<double> _x, std::vector<double> _y): nodes(_x.size()), x(_x), y(_y) {

	// Spline scrubbing: If two nodes are the same or not crescient, the node is erased
	scrubbing();

	// Create the gradients vector
	grad.resize(nodes - 1);

	for (int i = 0; i < nodes - 1; i++) {

		grad[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);

	}
	
	/// Initialize AlgLIB data
	alg_x.setcontent(x.size(), convert_data(x));
	alg_y.setcontent(x.size(), convert_data(y));

	/// Initialize limits
	x_min = _x[0];
	x_max = _x[_x.size() - 1];

	/// Create gradients vector


}



