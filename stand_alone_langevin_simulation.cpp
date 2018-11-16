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

#include "stand_alone_langevin_simulation.h"

StandAloneLangevinSimulation::StandAloneLangevinSimulation(raw_configuration_data _xml_data):
	StandAloneSimulation(_xml_data),
	V(_xml_data.volume),
	max_particles(_xml_data.max_part),
	min_particles(_xml_data.min_part),
	VTK_PATH(_xml_data.vtk_save),
	SAVE_VTK(_xml_data.SAVE_VTK_STEPS),
	max_dt(_xml_data.max_time_step) {

	// Print Simulation data
	std::cout << "Volume [m3]: " << V << std::endl;
	std::cout << "Max Particles [#]: " << max_particles << std::endl;
	std::cout << "Min Particles [#]: " << min_particles << std::endl;
	std::cout << "VTK Files Saving Path: " << VTK_PATH << std::endl;
	std::cout << "VTK Files Saving Frequency " << SAVE_VTK << std::endl;
	std::cout << "max dt allowed [sec]: " << max_dt << std::endl;

}

void StandAloneLangevinSimulation::run_simulation() {

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
									   
	// Simulation starting paramenters
	double T_start = this->start_T;
	double T_end = this->end_T;
	double time_end = this->end_Time;
	double V_start = this->V;
	double dt = 1.0e-09;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, m_fracts);

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(gas_species[0]);

	// setup the particle phase
	ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> pp(V_start);
	// set max number of particles
	pp.set_max_agg_number(this->max_particles);
	// set min num of partice
	pp.set_min_agg_number(this->min_particles);

	double elapsed_time = 0.0; // Simulation time
	int iterations = 0; // Iterations
	double snap_count = 0.0; // Counter for saving VTK file

	double T = T_start;

	WallClock clock;
	double SAVE_SNAPSHOT = this->SAVE_VTK;
	const int PRINT_EVERY =this->PRINT_STEPS;
	const int SAVE_EVERY = this->SAVE_STEPS;

	// paths where to save data
	std::ofstream plot_data;
	std::string plot_dir = this->SAVE_PATH;
	std::string filename = plot_dir + "Langevin_plot.dat";
	std::string vtk_path = this->VTK_PATH;

	// Plot data file open
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

	// PSD DATA
	std::ofstream part_sizes_file;
	std::string part_size_filename = plot_dir + "particles_sizes.dat";
	part_sizes_file.open(part_size_filename);
	// Aggregates diameters
	std::ofstream agg_sizes_file;
	std::string agg_size_filename = plot_dir + "aggregates_sizes.dat";
	agg_sizes_file.open(agg_size_filename);

	

	// Simulation Main Cycle
	while (elapsed_time < time_end) {

		// check the smallest particle in the system
		double d_min = pp.get_particles_smallest_diameter();
		// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
		double dt_max_lang = d_min*gas_species[0].get_bulk_density(gp.get_T()) / gp.get_gas_flux();
		// calculate dt max to have v*dt < d/2 for the smallest particle
		double dt_max_coll = sqrt(M_PI*gas_species[0].get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

		if (pp.get_aggregates_number() >= 1) {
			dt = std::min(dt_max_coll, dt_max_lang);
			if (dt >this->max_dt)
				dt = 1.0e-9;
			if (dt < 1.0e-12) {
				std::cout << "Stability: " << dt_max_lang << " Collision: " << dt_max_coll << " Smallest Diameter: " << d_min << std::endl;
			}
		}

		// particle phase timestep
		double g_si = pp.timestep_lc(dt, gp, cnt, T);

		// gas phase time step
		if (T <= this->end_T) dTdt = 0.0;
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		// Update elapsed time and iterations
		elapsed_time += dt;
		snap_count += dt;
		iterations++;

		if (counter_trigger(iterations, PRINT_EVERY)) {

			clock.stop();

			double T = gp.get_T();
			double S = gp.get_S("Si");

			std::cout
				<< "time: " << elapsed_time << '\t'									// time
				<< "Temp [K]: " << gp.get_T() << '\t'							// temperature
				<< "SuperSat: " << gp.get_S("Si") << '\t'						// supersaturation (S)
				<< "Nucl_Rate: " << cnt.nucleation_rate(T, S) << '\t'			// J
				<< "Gas-Phase: " << gp.get_n("Si") << '\t'						// ns
				<< "Stable_C_size: " << cnt.stable_cluster_size(T, S) << '\t'		// j
				<< "AVG_Part_per_Agg: " << pp.get_mean_particles_number() << '\t'		// N_m
				<< "AVG_Sint_level: " << pp.get_mean_sintering_level() << '\t'		//
				<< "AVG Sperical_Diam: " << pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< "Agg_Number: " << pp.get_aggregates_number() << '\t'
				<< "Agg_density: " << pp.get_aggregates_density() << '\t'
				<< "ts: " << dt << '\t'
				<< "Volume: " << pp.get_volume() << '\t'
				<< "Fractal_Dim: " << pp.get_mean_fractal_dimension() << '\t'
				<< "ts_exec_time: " << clock.interval() / PRINT_EVERY << std::endl;

			clock.start();
		}

		// Print file for plotting
		if (counter_trigger(iterations, SAVE_EVERY)) {

			clock.stop();

			double T = gp.get_T();
			double S = gp.get_S("Si");

			plot_data
				<< elapsed_time << '\t'
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

			clock.start();
		}

		// Save VTK
		if (snap_count >= SAVE_SNAPSHOT && pp.get_aggregates_number() > 0) {
			pp.save_vtk(iterations, vtk_path);
			snap_count = 0.0;
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

	system("PAUSE");

	return;

}
