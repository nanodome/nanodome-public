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

#ifndef SPECIES_H
#define SPECIES_H

#include "nanodome.h"
#include <string>
#include <cmath>


class Species {

// fundamental properties
    std::string formula; ///< Species formula (e.g. "ZnO")
    std::string name;    ///< IUPAC species name (e.g. "zinc monoxide")

    double mass; ///< Species mass [kg]

    double T_melt; ///< Melting temperature [K]

    double sigma; ///< L-J sigma value [m]
    double eps;   ///< L-J epsilon value [J]

// CNT oriented properties

    // densities
    double bulk_density_liq; ///< Density of the liquid bulk phase [kg/m3]
    double bulk_density_sol; ///< Density of the solid bulk phase [kg/m3]

    // surface tension is expressed as s_ten_A - s_ten_B*(T - s_ten_C) [N/m]
    double s_ten_A; ///< Surface tension coefficient A [N/m]
    double s_ten_B; ///< Surface tension coefficient B [N/(m K)]
    double s_ten_C; ///< Surface tension coefficient C [K]

    // saturation pressure is expressed as log10(psat) = (psat_A - psat_B/T) * 1.01e5 [Pa]
    double p_sat_A; ///< Saturation pressure coefficient A
    double p_sat_B; ///< Saturation pressure coefficient B [K]

public:

    /// Constructor
    /// \param _formula IUPAC name of the molecule
    Species(std::string _formula);

    /// Get IUPAC formula
    std::string get_formula() const { return formula; }

    /// Species mass [kg]
    double get_mass() const { return mass; }

    /// Species L-J sigma [m]
    double get_sigma() const { return sigma; }

    /// Species L-J epsilon [J]
    double get_epsilon() const { return eps; }

    /// Bulk material surface tension [N/m]
    /// \param T temperature [K]
	/*double s_ten(double T) const { 
		double s_ten_v = s_ten_A - s_ten_B*(T - s_ten_C); 
		return s_ten_v; }*/
	double s_ten(double T) const {
		double s_ten_v = (948.0 - 0.202*T)*0.001;
		return s_ten_v;
	} // OMLY FOR ALUMINUM

    /// Saturation pressure [Pa]
    /// \param T temperature [K]
    double p_sat(double T) const { return 1.01e5 * pow(10.0,(p_sat_A-(p_sat_B/T))); }
	//double p_sat(double T) const { return std::exp( 13.07-(36373.0/T))*1.01e5; } // ONLY FOR ALUMINUM

    /// Saturation density [#/m3]
    /// \param T temperature [K]
    double n_sat(double T) const { return p_sat(T)/(K_BOL*T); }

    /// Molecular volume based on liquid density [m3]
    double m_volume() const { return mass/(AMU*N_AVO*1000*bulk_density_liq); }

    /// Molecular surface based on liquid density and spherical assumption [m2]
    double m_surface() const { return 4.8360*pow(m_volume(),2./3.); }

    /// Get the bulk density depending on the temperature
    double get_bulk_density(double T) const { return (T<T_melt) ? bulk_density_sol : bulk_density_liq; }
};

#endif // SPECIES_H
