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

#ifndef STREAMLINE
#define STREAMLINE

#include <list>
#include <string>
#include <vector>
#include <iostream>

#include "species.h"
#include "splinelinear.h"

class streamline{

    /// Streamline ID
    int ID;

    /// Temperatures of the streamline in function of time [K]
    std::vector<double> T;

    /// Pressures of the streamline in funcion of time [Pa]
    std::vector<double> P;

    /// Time samples of the streamline [sec.]
    std::vector<double> Time;

    /// Species in the streamline in funtion of time
    std::list<Species> species;

    /// Species molar concentration for each time sample
    std::vector<double> Molar_Conc;

	/// Species molar concentration
	std::vector<double> Molar_Frac;

	/// Streamline coordinates
	std::vector<double> X;

	std::vector<double> Y;

	std::vector<double> Z;

	/// Linear Spline interpolating the Temperature samples
	SplineLinear* sT;

	/// Linear Spline interpolating the Pressure samples
	SplineLinear* sP;

	/// Linear Spline interpolating the Molar Concentration samples
	SplineLinear* sMC;

	/// Linear Spline interpolating the Molar Fraction samples
	SplineLinear* sMF;

	/// Linear Splines interpolating the Positions
	SplineLinear* sX; SplineLinear* sY; SplineLinear* sZ;

public:

    /// Default Constructor
    streamline();

    /// Parametric Constructor
    streamline(int _ID, std::vector<double> _T, std::vector<double> _P,
               std::vector<double> _Time, std::list<Species> _species,
			   std::vector<double> _Molar_Conc, std::vector<double> _Molar_Frac, 
			   std::vector<double> _X, std::vector<double> _Y, std::vector<double> _Z);

    /// Return Temperature [K] array
    std::vector<double> get_Temp() const {return T;}

    /// Return Pressure [Pa] array
    std::vector<double> get_Press() const {return P;}

    /// Return Times [sec]
    std::vector<double> get_Time() {return Time;}

	/// Return Molar Concentration for each species [sec]
	std::vector<double> get_Molar() const { return Molar_Conc; }

	/// Return Molar Fraction for each species [sec]
	std::vector<double> get_Molar_Fraction() const { return Molar_Frac; }

    /// Return Species
    std::list<Species> get_Species() const {return species;}

	/// Return Streamline ID
	int get_ID() const { return ID; }

	/// Return X Positions
	std::vector<double> get_x() const { return X; }

	/// Return Y Positions
	std::vector<double> get_y() const { return Y; }

	/// Return Z Positions
	std::vector<double> get_z() const { return Z; }

	/// Returns the Linear Spline interpolating the Temperature
	SplineLinear* get_splineT() { return sT; };

	/// Returns the Linear Spline interpolating the Pressure
	SplineLinear* get_splineP() { return sP; };

	/// Returns the Linear Spline interpolating the Molar Concentration
	SplineLinear* get_splineMC() { return sMC; };

	/// Returns the Linear Spline interpolating the Molar Fraction
	SplineLinear* get_splineMF() { return sMF; };

	/// Returns the Linear Spline interpolating the X coordinate
	SplineLinear* get_splineX() { return sX; };

	/// Returns the Linear Spline interpolating the X coordinate
	SplineLinear* get_splineY() { return sY; };

	/// Returns the Linear Spline interpolating the X coordinate
	SplineLinear* get_splineZ() { return sZ; };

	/// Returns the min time sample
	double get_minTime() { return Time[0]; };

	/// Returns the max time sample
	double get_maxTime() { return Time[Time.size() - 1]; };

    /// Print Streamline
    void printstream();

	/// Creates splines with data
	void create_splines();

};

#endif // STREAMLINE

