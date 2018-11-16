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

#include "species.h"

#include <iostream>
#include <cstdlib>


// for now only Si and Ar are implemented
Species::Species(std::string _formula) {

    if (_formula == "Si") {

        formula = "Si";
        name = "silicon";

        mass = 28.085*AMU;

        T_melt = 1687.;

        sigma = 3.30e-10;
        eps = 4.37e-20;

        bulk_density_liq = 2570.;
        bulk_density_sol = 2329.;

        s_ten_A = 0.732;
        s_ten_B = 0.000086;
        s_ten_C = 1685.;

        // validity range 300 - 3500 K
        p_sat_A = 7.5341;
        p_sat_B = 23399.;

    } else if (_formula == "Ar") {

        formula = "Ar";
        name = "argon";

        mass = 39.948*AMU;

        T_melt = 83.81;

        sigma = 3.41e-10;
        eps = 1.65e-21;

        bulk_density_liq = 1395.4;
        bulk_density_sol = 1395.4;

        s_ten_A = 0.;
        s_ten_B = 0.;
        s_ten_C = 0.;

        p_sat_A = 0.;
        p_sat_B = 0.;

    }
	else if (_formula == "Ti") {

		formula = "Ti";
		name = "titanium";

		mass = 47.867*AMU;

		T_melt = 1941.;

		sigma = 0.;
		eps = 0.;

		bulk_density_liq = 4110.0;
		bulk_density_sol = 4507.0;

		//s_ten_A = 1.557*1.0e+3;
		//s_ten_B = 0.156;
		//s_ten_C = 1941; //Ti menting point in K
		s_ten_A = 1.557;
		s_ten_B = 0.000156;
		s_ten_C = 1941; //Ti menting point in K

		// validity 300 - 4000
		p_sat_A = 15.5;
		p_sat_B = 55200;

	}
	else if (_formula == "O") {

		formula = "O";
		name = "oxygen";

		mass = 15.999*AMU;

		T_melt = 58.8;

		sigma = 0.;
		eps = 0.;

		bulk_density_liq = 1146.0;
		bulk_density_sol = 1146.0;

		s_ten_A = 0.;
		s_ten_B = 0.;
		s_ten_C = 0.;

		p_sat_A = 0.;
		p_sat_B = 0.;
	}
	else if (_formula == "Al") {
		formula = "Al";
		name = "aluminum";

		mass = 26.981*AMU;

		T_melt = 933.47;

		sigma = 0.;
		eps = 0.;

		bulk_density_liq = 2700.0;
		bulk_density_sol = 2300.0;

		s_ten_A = 0.;
		s_ten_B = 0.;
		s_ten_C = 0.;

		p_sat_A = 5.911;
		p_sat_B = 16211.0;
	}
	else {
		std::cout << "Species " << _formula << " does not exist!!!" << std::endl;
		std::exit(1);
	}
}
