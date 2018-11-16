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

#include "stand_alone_pbm_simulation.h"

StandAlonePBMSimulation::StandAlonePBMSimulation(raw_configuration_data _xml_data) :
	StandAloneSimulation(_xml_data),
	V(_xml_data.volume),
	max_particles(_xml_data.max_part),
	min_particles(_xml_data.min_part),
	fractal_dimension(_xml_data.frac_dim) {

	// Print Simulation data
	std::cout << "Volume [m3]: " << V << std::endl;
	std::cout << "Max Particles [#]: " << max_particles << std::endl;
	std::cout << "Min Particles [#]: " << min_particles << std::endl;
	std::cout << "Fractal Dimension: " << fractal_dimension << std::endl;


}

void StandAlonePBMSimulation::run_simulation() {

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

	// gas phase info
	double p = this->pressure; // pressure [Pa]
	double dTdt = this->temp_gradient; // temperature gradient [K/s]

	double T_start = this->start_T;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, m_fracts);

	// setup the nucleation theory we want to use (for the moment monospecies)
	ClassicalNucleationTheory cnt(gas_species[0]);

	// setup the particle phase we want to use
	PBMFractalParticlePhase<PBMAggregate<Particle>> pp(this->fractal_dimension, this->V);

	// set max and min number of aggregates
	pp.set_max_aggregates(this->max_particles);
	pp.set_min_aggregates(this->min_particles);

	double t = 0.0;
	int iter = 0;


	int PRINT_EVERY = this->PRINT_STEPS;

	int SAVE_EVERY = this->SAVE_STEPS;

	// Plot Data
	std::ofstream plot_data;
	std::string dir = this->SAVE_PATH;
	std::string filename = dir + "PBM_plot.dat";
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

	// PSD Data
	// Prinary particles
	std::ofstream part_sizes_file;
	std::string part_size_filename = dir + "particles_sizes.dat";
	part_sizes_file.open(part_size_filename);
	// Aggregates diameters
	std::ofstream agg_sizes_file;
	std::string agg_size_filename = dir + "aggregates_sizes.dat";
	agg_sizes_file.open(agg_size_filename);

	double END_TIME = this->end_Time;

	clock.start();

	// loop over timesteps
	while (t < END_TIME) {

		if (gp.get_T() < this->end_T){
			dTdt = 0;}

		// species source term for the gas phase
		double g_si = 0.0;

		// calculate the timestep using an exponential waiting time
		double R_tot = pp.get_total_processes_rate(gp, cnt);
		double rho = ndm::uniform_double_distr(ndm::rand_gen);

		// exponential waiting time
		double dt = -log(rho) / R_tot;

		// Strang first step
		gp.timestep(dt / 2.0, dTdt, 0);
		pp.volume_expansion(dt / 2.0, gp);

		// Strang second step
		g_si += pp.timestep(dt, gp, cnt);
		gp.timestep(dt, 0, 0, { -g_si,0.0 });

		// Strang third step
		gp.timestep(dt / 2.0, dTdt, 0);
		pp.volume_expansion(dt / 2.0, gp);

		t += dt;
		iter++;

		if (counter_trigger(iter, PRINT_EVERY)) {

			clock.stop();

			double T = gp.get_T();
			double S = gp.get_S("Si");

			std::cout 
				<< "time: " << t << '\t'									// time
				<< "Temp [K]: " << gp.get_T() << '\t'							// temperature
				<< "SuperSat: " << gp.get_S("Si") << '\t'						// supersaturation (S)
				<< "Nucl_Rate: " << cnt.nucleation_rate(T, S) << '\t'			// J
				<< "Gas-Phase: " << gp.get_n("Si") << '\t'						// ns
				<< "Stable_C_size: " << cnt.stable_cluster_size(T, S) << '\t'		// j
				<< "AVG_Part_per_Agg[m]: " << pp.get_mean_particles_number() << '\t'		// N_m
				<< "AVG_Sint_level[%]: " << pp.get_mean_sintering_level() << '\t'		//
				<< "AVG Sperical_Diam[m]: " << pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< "Agg_Number[#]: " << pp.get_aggregates_number() << '\t'
				<< "Agg_density[#/m3]: " << pp.get_aggregates_density() << '\t'
				<< "Volume[m3]: " << pp.get_volume() << '\t'
				<< "Fractal_Dim: " << pp.get_mean_fractal_dimension() << '\t'
				<< "ts_exec_time[sec]: " << clock.interval() / PRINT_EVERY << std::endl;

			clock.start();
		}

		// Print file for plotting
		if (counter_trigger(iter, SAVE_EVERY)) {

			clock.stop();

			double T = gp.get_T();
			double S = gp.get_S("Si");

			plot_data
				<< t << '\t'
				<< gp.get_T() << '\t'
				<< gp.get_S("Si") << '\t'
				<< cnt.nucleation_rate(T, S) << '\t'
				<< gp.get_n("Si") << '\t'
				<< cnt.stable_cluster_size(T, S) << '\t'
				<< pp.get_mean_particles_number() << '\t'
				<< pp.get_mean_sintering_level() << '\t'
				<< pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< pp.get_aggregates_number() << '\t'
				<< pp.get_aggregates_density() << '\t'
				<< pp.get_volume() << '\t'
				<< pp.get_mean_fractal_dimension() << '\t'
				<< clock.interval() / PRINT_EVERY << std::endl;

			// Print file

			clock.start();
		}
	}
	// Close plot file
	plot_data.close();


	// Anyway at the end save aggregates data for PSD
	std::valarray<double> particles_sizes = pp.get_particles_sizes();
	std::valarray<double> aggregates_sizes = pp.get_aggregates_sizes();

	// print particles sizes
	for (int i = 0; i < particles_sizes.size(); i++) {
		part_sizes_file << particles_sizes[i] << " ";
	}
	part_sizes_file << std::endl;
	// print aggregates sizes
	for (int i = 0; i < aggregates_sizes.size(); i++) {
		agg_sizes_file << aggregates_sizes[i] << " ";
	}
	agg_sizes_file << std::endl;

	part_sizes_file.close();
	agg_sizes_file.close();


	return;

}
