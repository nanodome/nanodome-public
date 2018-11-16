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

#include "xml_templates.h"

void linked_moments_template(std::string _filename) {

	std::ofstream xml_out;
	xml_out.open(_filename);

	xml_out << "<MESO_SIM>" << std::endl;
	xml_out << "<CONDITIONS>" << std::endl;
	xml_out << "<LINK_PATH temp_grad = yes / no(attribute for using only the temperature gradient or all the streamline data)>path to NanoDome XML streamlines file</LINK_PATH>" << std::endl;
	xml_out << "</CONDITIONS>" << std::endl;

	xml_out << "<DATA_PLOT>" << std::endl;
	xml_out << "<SAVE_STEPS>steps frequency for saving data to plot</SAVE_STEPS>" << std::endl;
	xml_out << "<PRINT_STEPS>steps frequency for printing on screen data</PRINT_STEPS>" << std::endl;
	xml_out << "<SAVE_PATH>Path for saving data</SAVE_PATH>" << std::endl;
	xml_out << "</DATA_PLOT>" << std::endl;

	xml_out << "<MODEL>" << std::endl;
	xml_out << "<MODEL_TYPE>Moments</MODEL_TYPE>" << std::endl;
	xml_out << "<TIME_STEP>Simulation Timestep</TIME_STEP>" << std::endl;
	xml_out << "</MODEL>" << std::endl;
	xml_out << "</MESO_SIM>" << std::endl;

	std::cout << "TEMPLATE CREATED!!" << std::endl;

	xml_out.close();
	return;
}
void linked_pbm_template(std::string _filename) {

	std::ofstream xml_out;
	xml_out.open(_filename);

	xml_out << "<MESO_SIM>" << std::endl;
	xml_out << "<CONDITIONS>" << std::endl;
	xml_out << "<LINK_PATH temp_grad = yes / no(attribute for using only the temperature gradient or all the streamline data)>path to NanoDome XML streamlines file</LINK_PATH>" << std::endl;
	xml_out << "</CONDITIONS>" << std::endl;

	xml_out << "<DATA_PLOT>" << std::endl;
	xml_out << "<SAVE_STEPS>steps frequency for saving data to plot</SAVE_STEPS>" << std::endl;
	xml_out << "<PRINT_STEPS>steps frequency for printing on screen data</PRINT_STEPS>" << std::endl;
	xml_out << "<SAVE_PATH>Path for saving data</SAVE_PATH>" << std::endl;
	xml_out << "</DATA_PLOT>" << std::endl;

	xml_out << "<MODEL>" << std::endl;
	xml_out << "<MODEL_TYPE>PBM</MODEL_TYPE>" << std::endl;
	xml_out << "<VOLUME>Control Volume in m3</VOLUME>" << std::endl;
	xml_out << "<MAX_PART>Maximum number of particles</MAX_PART>" << std::endl;
	xml_out << "<MIN_PART>Minimim number of particles</MIN_PART>" << std::endl;
	xml_out << "<F_DIM>Nanoparticles fractal dimension</F_DIM>" << std::endl;
	xml_out << "</MODEL>" << std::endl;
	xml_out << "</MESO_SIM>" << std::endl;

	std::cout << "TEMPLATE CREATED!!" << std::endl;

	xml_out.close();
	return;
}

void linked_langevin_template(std::string _filename) {

	std::ofstream xml_out;
	xml_out.open(_filename);

	xml_out << "<MESO_SIM>" << std::endl;
	xml_out << "<CONDITIONS>" << std::endl;
	xml_out << "<LINK_PATH temp_grad = yes / no(attribute for using only the temperature gradient or all the streamline data)>path to NanoDome XML streamlines file</LINK_PATH>" << std::endl;
	xml_out << "</CONDITIONS>" << std::endl;

	xml_out << "<DATA_PLOT>" << std::endl;
	xml_out << "<SAVE_STEPS>steps frequency for saving data to plot</SAVE_STEPS>" << std::endl;
	xml_out << "<PRINT_STEPS>steps frequency for printing on screen data</PRINT_STEPS>" << std::endl;
	xml_out << "<SAVE_VTK_STEPS>time frequency for saving data for visualization (VTK)</SAVE_VTK_STEPS>" << std::endl;
	xml_out << "<SAVE_VTK_PATH>Path for saving visualization data</SAVE_VTK_PATH>" << std::endl;
	xml_out << "<SAVE_PATH>Path for saving data</SAVE_PATH>" << std::endl;
	xml_out << "</DATA_PLOT>" << std::endl;

	xml_out << "<MODEL>" << std::endl;
	xml_out << "<MODEL_TYPE>Langevin</MODEL_TYPE>" << std::endl;
	xml_out << "<VOLUME>Control Volume in m3</VOLUME>" << std::endl;
	xml_out << "<MAX_PART>Maximum number of particles</MAX_PART>" << std::endl;
	xml_out << "<MIN_PART>Minimim number of particles</MIN_PART>" << std::endl;
	xml_out << "<MAX_TS>Maximum timestep feaseble for the method stability</MAX_TS>"<<std::endl;
	xml_out << "</MODEL>" << std::endl;

	xml_out << "</MESO_SIM>" << std::endl;

	std::cout << "TEMPLATE CREATED!!" << std::endl;

	xml_out.close();
	return;
}

void stand_alone_moment_template(std::string _filename) {

	std::ofstream xml_out;
	xml_out.open(_filename);

	xml_out << "<MESO_SIM>" << std::endl;
	xml_out << "<CONDITIONS>" << std::endl;
	xml_out << "<C_SPECIES>Condensating species</C_SPECIES>" << std::endl;
	xml_out << "<MOLAR_F_C>Condensating species starting molar fraction int %< / MOLAR_F_C>" << std::endl;
		xml_out << "<B_SPECIES>Bath species</B_SPECIES>" << std::endl;
		xml_out << "<MOLAR_F_B>Bath species molar fraction in %</MOLAR_F_B>" << std::endl;
		xml_out << "<PRESSURE>pressure in Pa< / PRESSURE>" << std::endl;
		xml_out << "<S_TEMP>Simulation starting temperature in K</S_TEMP>" << std::endl;
		xml_out << "<E_TEMP>Simulation ending temperature in K</E_TEMP>" << std::endl;
		xml_out << "<END_TIME>Simulation ending time in sec</END_TIME>" << std::endl;
		xml_out << "<T_GRAD>Temperature gradient</T_GRAD>" << std::endl;
	xml_out << "</CONDITIONS>" << std::endl;

		xml_out << "<DATA_PLOT>" << std::endl;
		xml_out << "<SAVE_STEPS>steps frequency for saving data to plot</SAVE_STEPS>" << std::endl;
		xml_out << "<PRINT_STEPS>steps frequency for printing on screen data</PRINT_STEPS>" << std::endl;
		xml_out << "<SAVE_PATH>Path for saving data</SAVE_PATH>" << std::endl;
		xml_out << "</DATA_PLOT>" << std::endl;

		xml_out << "<MODEL>" << std::endl;
		xml_out << "<MODEL_TYPE>Moments</MODEL_TYPE>" << std::endl;
		xml_out << "<TIME_STEP>Simulation Timestep</TIME_STEP>" << std::endl;
		xml_out << "</MODEL>" << std::endl;

		xml_out << "</MESO_SIM>" << std::endl;

		std::cout << "TEMPLATE CREATED!!" << std::endl;

	xml_out.close();
	return;
}

void stand_alone_pbm_template(std::string _filename) {

	std::ofstream xml_out;
	xml_out.open(_filename);

	xml_out << "<MESO_SIM>" << std::endl;
	xml_out << "<CONDITIONS>" << std::endl;
	xml_out << "<C_SPECIES>Condensating species</C_SPECIES>" << std::endl;
	xml_out << "<MOLAR_F_C>Condensating species starting molar fraction int %< / MOLAR_F_C>" << std::endl;
	xml_out << "<B_SPECIES>Bath species</B_SPECIES>" << std::endl;
	xml_out << "<MOLAR_F_B>Bath species molar fraction in %</MOLAR_F_B>" << std::endl;
	xml_out << "<PRESSURE>pressure in Pa< / PRESSURE>" << std::endl;
	xml_out << "<S_TEMP>Simulation starting temperature in K</S_TEMP>" << std::endl;
	xml_out << "<E_TEMP>Simulation ending temperature in K</E_TEMP>" << std::endl;
	xml_out << "<END_TIME>Simulation ending time in sec</END_TIME>" << std::endl;
	xml_out << "<T_GRAD>Temperature gradient</T_GRAD>" << std::endl;
	xml_out << "</CONDITIONS>" << std::endl;

	xml_out << "<DATA_PLOT>" << std::endl;
	xml_out << "<SAVE_STEPS>steps frequency for saving data to plot</SAVE_STEPS>" << std::endl;
	xml_out << "<PRINT_STEPS>steps frequency for printing on screen data</PRINT_STEPS>" << std::endl;
	xml_out << "<SAVE_PATH>Path for saving data</SAVE_PATH>" << std::endl;
	xml_out << "</DATA_PLOT>" << std::endl;

	xml_out << "<MODEL>" << std::endl;
	xml_out << "<MODEL_TYPE>PBM</MODEL_TYPE>" << std::endl;
	xml_out << "<VOLUME>Control Volume in m3</VOLUME>" << std::endl;
	xml_out << "<MAX_PART>Maximum number of particles</MAX_PART>" << std::endl;
	xml_out << "<MIN_PART>Minimim number of particles</MIN_PART>" << std::endl;
	xml_out << "<F_DIM>Nanoparticles fractal dimension</F_DIM>" << std::endl;
	xml_out << "</MODEL>" << std::endl;
	xml_out << "</MESO_SIM>" << std::endl;

	std::cout << "TEMPLATE CREATED!!" << std::endl;

	xml_out.close();
	return;
}

void stand_alone_langevin_template(std::string _filename) {

	std::ofstream xml_out;
	xml_out.open(_filename);

	xml_out << "<MESO_SIM>" << std::endl;
	xml_out << "<CONDITIONS>" << std::endl;
	xml_out << "<C_SPECIES>Condensating species</C_SPECIES>" << std::endl;
	xml_out << "<MOLAR_F_C>Condensating species starting molar fraction int %< / MOLAR_F_C>" << std::endl;
	xml_out << "<B_SPECIES>Bath species</B_SPECIES>" << std::endl;
	xml_out << "<MOLAR_F_B>Bath species molar fraction in %</MOLAR_F_B>" << std::endl;
	xml_out << "<PRESSURE>pressure in Pa< / PRESSURE>" << std::endl;
	xml_out << "<S_TEMP>Simulation starting temperature in K</S_TEMP>" << std::endl;
	xml_out << "<E_TEMP>Simulation ending temperature in K</E_TEMP>" << std::endl;
	xml_out << "<END_TIME>Simulation ending time in sec</END_TIME>" << std::endl;
	xml_out << "<T_GRAD>Temperature gradient</T_GRAD>" << std::endl;
	xml_out << "</CONDITIONS>" << std::endl;

	xml_out << "<DATA_PLOT>" << std::endl;
	xml_out << "<SAVE_STEPS>steps frequency for saving data to plot</SAVE_STEPS>" << std::endl;
	xml_out << "<PRINT_STEPS>steps frequency for printing on screen data</PRINT_STEPS>" << std::endl;
	xml_out << "<SAVE_VTK_STEPS>time frequency for saving data for visualization (VTK)</SAVE_VTK_STEPS>" << std::endl;
	xml_out << "<SAVE_VTK_PATH>Path for saving visualization data</SAVE_VTK_PATH>" << std::endl;
	xml_out << "<SAVE_PATH>Path for saving data</SAVE_PATH>" << std::endl;
	xml_out << "</DATA_PLOT>" << std::endl;

	xml_out << "<MODEL>" << std::endl;
	xml_out << "<MODEL_TYPE>Langevin</MODEL_TYPE>" << std::endl;
	xml_out << "<VOLUME>Control Volume in m3</VOLUME>" << std::endl;
	xml_out << "<MAX_PART>Maximum number of particles</MAX_PART>" << std::endl;
	xml_out << "<MIN_PART>Minimim number of particles</MIN_PART>" << std::endl;
	xml_out << "<MAX_TS>Maximum timestep feaseble for the method stability</MAX_TS>" << std::endl;
	xml_out << "</MODEL>" << std::endl;

	xml_out << "</MESO_SIM>" << std::endl;

	std::cout << "TEMPLATE CREATED!!" << std::endl;

	xml_out.close();
	return;
}
