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

#include "streamline.h"

streamline::streamline(){}

streamline::streamline(int _ID, std::vector<double> _T, std::vector<double> _P,
	std::vector<double> _Time, std::list<Species> _species,
	std::vector<double> _Molar_Conc, std::vector<double> _Molar_Frac,
	std::vector<double> _X, std::vector<double> _Y, std::vector<double> _Z) :
	ID(_ID), T(_T), P(_P), Time(_Time),
	species(_species), Molar_Conc(_Molar_Conc), Molar_Frac(_Molar_Frac),
	X(_X), Y(_Y), Z(_Z)
{}

void streamline::create_splines() {

	// Create Splines

	// Temperature
	if (Time.size() > 1 && T.size() > 1)
		sT = new SplineLinear(Time, T);
	else {
		/// Error
	}

	// Pressure
	if (Time.size() > 1 && P.size() > 1)
		sP = new SplineLinear(Time, P);
	else{
		///error
	}

	// Molar Concentration
	if (Time.size() > 1 && Molar_Conc.size() > 1) {
		sMC = new SplineLinear(Time, Molar_Conc);
	}
	else{
		// error
	}

	// Molar Fraction
	if (Time.size() > 1 && Molar_Frac.size() > 1) {
		sMF = new SplineLinear(Time, Molar_Frac);
	}
	else{
		// Error
	}

	// X
	if (Time.size() > 1 && X.size() > 1) {
		sX = new SplineLinear(Time, X);
	}
	else
	{
		//error
	}

	// Y
	if (Time.size() > 1 && Y.size() > 1)
		sY = new SplineLinear(Time, Y);
	else{
		// Error
	}

	// Z 
	if (Time.size() > 1 && Z.size() > 1)
		sZ = new SplineLinear(Time, Z);
	else {
		// error
	}

}

    
void streamline::printstream() {

	std::cout << "ID: " << ID << std::endl;

	std::cout << "Time Samples" << std::endl;
	for (auto t = Time.begin(); t != Time.end(); ++t) std::cout << (*t) << " ";
	std::cout << std::endl;

	std::cout << "Temperatures" << std::endl;
	for (auto t = T.begin(); t != T.end(); ++t) std::cout << (*t) << " ";
	std::cout << std::endl;

	std::cout << "Pressure" << std::endl;
	for (auto t = P.begin(); t != P.end(); ++t) std::cout << (*t) << " ";
	std::cout << std::endl;

	std::cout << "Species" << std::endl;
	for (auto t = species.begin(); t != species.end(); ++t)
		std::cout << (*t).get_formula() << " ";
	std::cout << std::endl;

	std::cout << "Molar Concentration" << std::endl;
	for (auto t = Molar_Conc.begin(); t != Molar_Conc.end(); ++t) std::cout << (*t) << " ";
	std::cout << std::endl;

	std::cout << "Molar Fraction" << std::endl;
	for (auto t = Molar_Frac.begin(); t != Molar_Frac.end(); ++t) std::cout << (*t) << " ";
	std::cout << std::endl;


	std::cout << std::endl;


}
