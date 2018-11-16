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

#include "stand_alone_moments_simulation.h"

StandAloneMomentsSimulation::StandAloneMomentsSimulation(raw_configuration_data _xml_data) :
	StandAloneSimulation(_xml_data), time_step(_xml_data.dt) {

	// Print Simulation data
	std::cout << "Moments Timestep: " << time_step << std::endl;


}

void StandAloneMomentsSimulation::run_simulation(){

	WallClock clock;

	// create the vector of species for the gas phase
	std::vector<Species> gas_species;

	// Condensing Species
	for (auto it = this->c_species.begin(); it != this->c_species.end(); it++) {
		Species s = Species((*it));
		gas_species.push_back(s);
	}

	// Carrier gas Species
	for (auto it = this->b_species.begin(); it != this->b_species.end(); it++) {
		Species s = Species((*it));
		gas_species.push_back(s);
	}
	// create the vector of molar fractions for each species (ordered following the species order)
	std::valarray<double> m_fracts(gas_species.size());
	int idx = 0;
	for (auto it = this->c_species_m_fraction.begin(); it != this->c_species_m_fraction.end(); it++) {
		m_fracts[idx] = (*it);
		idx++;
	}

	for (auto it = this->b_species_m_fraction.begin(); it != this->b_species_m_fraction.end(); it++) {
		m_fracts[idx] = (*it);
		idx++;
	}

	// simulation timestep
	const double dt = this->time_step;

	// gas phase info
	double p = this->pressure; // pressure [Pa]
	double dTdt = this->temp_gradient; // temperature gradient [K/s]

	double T_start = this->start_T;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, m_fracts);

	// setup the nucleation theory we want to use (for the moment monospecies)
	ClassicalNucleationTheory cnt(gas_species[0]);

	// set the moment method with the condensing species (for the moment monospecies)
	//MomentModelFriedlander mm(si);
	MomentModelPratsinis mm(gas_species[0]);

	double t = 0.;
	int iter = 0;

	int PRINT_EVERY = this->PRINT_STEPS;

	int SAVE_EVERY = this->SAVE_STEPS;

	std::ofstream plot_data;
	std::string dir = this->SAVE_PATH;
	std::string filename = SAVE_PATH + "MOMENTS_plot.dat";
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

	double END_TIME = this->end_Time;

	// loop over timesteps until the final temperature
	while (t < END_TIME) {

		if (gp.get_T()<this->end_T) dTdt = 0;
		
		double T = gp.get_T();
		double ns = gp.get_n("Si");
		double S = ns / gas_species[0].n_sat(T); // monospecies -> first element of the list
		double J = cnt.nucleation_rate(T, S);
		double j = cnt.stable_cluster_size(T, S);

		// moment method timestep and species consumption retrieval
		double g_si = mm.timestep(dt, gp, cnt);

		// updating the gas phase
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		++iter;
		t += dt;

		if (counter_trigger(iter, PRINT_EVERY)) {
			//gp.print();

			std::cout
				<< "time: " << t << '\t'
				<< "Temp [K]: " << T << '\t'
				<< "SuperSat: " << S << '\t'
				<< "Nucl_Rate [#/m3 sec]: " << J << '\t'
				<< "Gas-Phase [#]: " << ns << '\t'
				<< "Stable_C_size[#]: " << j << '\t'
				<< "Agg_density[#/m3]: " << mm.get_density() << '\t'
				<< "AVG Sperical_Diam[m]: " << mm.get_mean_diameter() << std::endl;
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
	return;

}
