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

#include "linked_moments_simulation.h"

LinkedMomentsSimulation::LinkedMomentsSimulation(raw_configuration_data _xml_data) :
	LinkedSimulation(_xml_data), time_step(_xml_data.dt) {

	// Print Simulation data
	std::cout << "Moments Timestep: " << time_step<<std::endl;

}

void LinkedMomentsSimulation::run_simulation() {

	WallClock clock;

	
	// Declare the path for the XML file
	//std::string xml_path = "streams/NanoDomeStreamlines.xml";
	std::string xml_path = this->STREAMS_PATH;

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	

	// Main cycle: a simulation for each streamline extracted from the CFD

	for (int s = 0; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;

		std::list<Species> l_gas_species = l_streams_xml[s].get_Species();
		std::vector<Species> gas_species;
		for (auto sp = l_gas_species.begin(); sp != l_gas_species.end(); sp++) {
			gas_species.push_back((*sp));
		}

		// Create splines interpolating the samples extracted
		l_streams_xml[s].create_splines();

		std::vector<Spline*> l_sc;
		l_sc.push_back(l_streams_xml[s].get_splineMF());

		// Dummy Spline for Ar
		SplineLinear *sC_ar = new SplineLinear();
		l_sc.push_back(sC_ar);

		// simulation timestep
		double dt = this->time_step;

		GasPhaseLink gp(l_streams_xml[s].get_splineT(),
			l_streams_xml[s].get_splineP(),
			gas_species, l_sc,
			l_streams_xml[s].get_minTime());


		// setup the nucleation theory we want to use
		ClassicalNucleationTheory cnt(gas_species[0]);

		// set the moment method with the condensing species
		//MomentModelFriedlander mm(si);
		MomentModelPratsinis mm(gas_species[0]);

		double t = l_streams_xml[s].get_minTime();

		// simulation end time
		double end_time = l_streams_xml[s].get_maxTime();

		int iter = 0;

		int PRINT_EVERY = this->PRINT_STEPS;

		int SAVE_EVERY = this->SAVE_STEPS;

		std::ofstream plot_data;
		
		std::string dir = this->SAVE_PATH;
		std::string filename = dir + "MOMENTS_plot" + std::to_string(s) + "Streamline.dat";
		plot_data.open(filename);

		// Print headlines
		plot_data 
			<< "Time[sec]" << '\t'
			<< "Temp[K]" << '\t'
			<< "SuperSaturation ratio" << '\t'
			<< "Nucl Rate" << '\t'
			<< "Species # density" << '\t'
			<< "Stable cluster size[m]" << '\t'
			<< "AVG Part Num[#]: " << '\t'
			<< "Sint level[%]" << '\t'
			<< "AVG diameter[m]" << '\t'
			<< "Agg. #[#]" << '\t'
			<< "Agg density[#/m3]" << '\t'
			<< "Volume[m]" << '\t'
			<< "AVG fract dim" << '\t'
			<< "ts exec time"
			<< std::endl;
		

		// loop over timesteps until the final temperature
		while (t < end_time) {

			double T = gp.get_T();
			double ns = gp.get_n("Si");
			double S = ns / gas_species[0].n_sat(T);
			double J = cnt.nucleation_rate(T, S);
			double j = cnt.stable_cluster_size(T, S);

			// moment method timestep and species consumption retrieval
			double g_si = mm.timestep(dt, gp, cnt);

			// updating the gas phase
			//gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });
			gp.timestep(t, dt);

			++iter;
			t += dt;

			// Step Print
			if (counter_trigger(iter, PRINT_EVERY)) {
				//gp.print();

				std::cout
					<< "time: "			<< t << '\t'
					<< "Temp: "				<< T << '\t'
					<< "SuperSat: "				<< S << '\t'
					<< "Nucl_Rate: " << J << '\t'
					<< "Gas-Phase: "		<< ns << '\t'
					<< "Stable_C_size: "		<< j << '\t'
					<< "Agg_density: "	<< mm.get_density() << '\t'
					<< "AVG Sperical_Diam: " << mm.get_mean_diameter() << std::endl;
			}

			// File savings
			if (counter_trigger(iter, SAVE_EVERY)) {
				plot_data
					<< t << '\t'
					<< T << '\t'
					<< S << '\t'
					<< J << '\t'
					<< ns << '\t'
					<< j << '\t'
					<< 0.0 << '\t' //dummy particles number
					<< 0.0 << '\t' // dummy sintering level
					<< mm.get_mean_diameter() << '\t'
					<< 0.0 << '\t' // dummy aggregates number
					<< mm.get_density() << '\t'
					<< 0.0 << '\t' // dummy Volume
					<< 0.0 << '\t' // dummy fractal dimension
					<< 0.0 << std::endl; // dummy interval time
			}
		}

		plot_data.close();

	}

	system("PAUSE");


	return;

}
