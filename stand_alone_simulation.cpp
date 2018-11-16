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

#include "stand_alone_simulation.h"

StandAloneSimulation::StandAloneSimulation(raw_configuration_data _xml_data) :
	MesoSimulation(_xml_data),
	c_species(_xml_data.c_species),
	b_species(_xml_data.b_species),
	c_species_m_fraction(_xml_data.c_s_molar_fraction),
	b_species_m_fraction(_xml_data.b_s_molar_fraction),
	pressure(_xml_data.pressure),
	start_T(_xml_data.start_temp),
	end_T(_xml_data.end_temp),
	end_Time(_xml_data.end_time),
	temp_gradient(_xml_data.t_grad) {

	// Print Data
	std::cout << "Nucleating Species: ";
	for (auto it : c_species) 
		std::cout << it;
	std::cout << std::endl;

	std::cout << "Bath Species: ";
	for (auto it : b_species)
		std::cout << it;
	std::cout << std::endl;

	std::cout << "Nucleating Species Molar Fractions: ";
	for (auto it : c_species_m_fraction)
		std::cout << it;
	std::cout << std::endl;

	std::cout << "Bath Species Molar Fractions: ";
	for (auto it : b_species_m_fraction)
		std::cout << it;
	std::cout << std::endl;

	std::cout << "Pressure [pa]: " << pressure << std::endl;

	std::cout << "Starting Temperature [K]: " << start_T << std::endl;

	std::cout << "Ending Temperature [K]: " << end_T << std::endl;

	std::cout << "Ending Time[sec]: " << end_Time << std::endl;

	std::cout << "Temperature Gradiend [K/sec]: " << temp_gradient << std::endl;



}
