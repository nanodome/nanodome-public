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

#ifndef _SPLINELINEAR_H
#define _SPLINELINEAR_H

#include "spline.h"

class SplineLinear : public Spline {

private:

    /// AlgLIB interpolator
    alglib::spline1dinterpolant s;

public:

	/// Blank Constructor
	SplineLinear();

    /// Constructor
    ///	\param std::vector<double> _x set of x values for the nodes
    /// \param std::vector<double> _y set of y values for the nodes
    SplineLinear(std::vector<double> _x, std::vector<double> _y);

    /// Get interpolated value
    double get_v(double _x);
};

#endif

