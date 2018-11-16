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

#ifndef _SPLINE_H
#define _SPLINE_H

// Including the interpolation library
#include "alglib/src/interpolation.h"

// STD
#include <vector>

// Error Management
#include "i_o.h"


class Spline {

protected:

    /// Minimum admissible value for x
    double x_min;

    /// Maximum admissible value for x
    double x_max;

    /// Number of Nodes
    int nodes;

    /// Set of nodes x
    std::vector<double> x;

    /// Set of node f(x) = y
    std::vector<double> y;

	/// Set of gradients among two nodes
	std::vector<double> grad;

    /// Library data arrays
    alglib::real_1d_array alg_x;
    alglib::real_1d_array alg_y;

	/// Function for scrubbing not coherent data in the x vector
	void scrubbing();

public:

    /// Blank Constructor
	Spline();

    /// Constructor
    ///	\param std::vector<double> _x set of x values for the nodes
    /// \param std::vector<double> _y set of y values for the nodes
    Spline(std::vector<double> _x, std::vector<double> _y);

    /// Get the interpolated value
    virtual double get_v(double _xi) = 0;

	/// Get the x minimum value
	double get_x_min() const { return x_min; }

	/// Get the x maximum value
	double get_x_max() const { return x_max; }

	/// Get x nodes
	std::vector<double> get_nodes_x()const { return x; }

	/// Get y nodes
	std::vector<double> get_nodes_y()const { return y; }

	/// Get gradients
	std::vector<double> get_gradients()const { return grad; }
};

#endif
