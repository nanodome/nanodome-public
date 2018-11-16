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

#include "linked_pbm_temp_simulation.h"

PBMLinkedTempSimulation::PBMLinkedTempSimulation(raw_configuration_data _xml_data) :
	LinkedSimulation(_xml_data),
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

void PBMLinkedTempSimulation::run_simulation() {

	WallClock clock;

	// Declare the path for the XML file
	std::string xml_path = this->STREAMS_PATH;

	// Global PSD from all streamlines (primary particles)
	std::ofstream global_psd_primary;
	std::string global_psd_primary_file = this->SAVE_PATH + "global_primary_PSD.dat";
	global_psd_primary.open(global_psd_primary_file);

	// Global PSD from all streamlines (aggregates)
	std::ofstream global_psd_aggregates;
	std::string global_psd_aggregates_file = this->SAVE_PATH + "global_aggregates_PSD.dat";
	global_psd_aggregates.open(global_psd_aggregates_file);

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

		// streamline start time
		double start_time = l_streams_xml[s].get_minTime();
		// starting molar fraction value
		double start_mf = l_streams_xml[s].get_splineMF()->get_v(start_time);

		GasPhaseLink gp(l_streams_xml[s].get_splineT(),
			l_streams_xml[s].get_splineP(),
			gas_species, start_time, start_mf);
		
		// setup the nucleation theory we want to use
		ClassicalNucleationTheory cnt(gas_species[0]);

		// setup the particle phase we want to use
		PBMFractalParticlePhase<PBMAggregate<Particle>> pp(this->fractal_dimension, this->V);

		// set max and min number of aggregates
		pp.set_max_aggregates(this->max_particles);
		pp.set_min_aggregates(this->min_particles);

		// get the streamline starting time
		double t = l_streams_xml[s].get_minTime();

		// get the end time for the streamline
		double end_time = l_streams_xml[s].get_maxTime();

		int iter = 0;

		const int PRINT_EVERY = this->PRINT_STEPS;

		const int SAVE_EVERY = this->SAVE_STEPS;
		std::ofstream plot_data;
		std::string dir = this->SAVE_PATH;
		std::string filename = dir + "PBM_plot" + std::to_string(s) + "Streamline.dat";
		plot_data.open(filename);

		// Print headlines
		plot_data << "Time[sec]" << '\t'
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

		// Prinary particles
		std::ofstream part_sizes_file;
		std::string part_size_filename = dir + "particles_sizes" + std::to_string(s) + "Streamline.dat";
		part_sizes_file.open(part_size_filename);
		// Aggregates diameters
		std::ofstream agg_sizes_file;
		std::string agg_size_filename = dir + "aggregates_sizes" + std::to_string(s) + "Streamline.dat";
		agg_sizes_file.open(agg_size_filename);


		clock.start();

		// loop over timesteps

		while (t < end_time) {

			// species source term for the gas phase
			double g_si = 0.0;

			// calculate the timestep using an exponential waiting time
			double R_tot = pp.get_total_processes_rate(gp, cnt);
			double rho = ndm::uniform_double_distr(ndm::rand_gen);
			// if RHO is equal to 0.0, is set to the minimum available number for double type
			if (rho == 0.0) {
				std::cout << "RHO = 0, -> Min double!" << std::endl;
				rho = std::numeric_limits<double>::min();
			}
			// if RHO is equal to 1 is subtracted an EPSILON value
			if (rho == 1.0) {
				std::cout << "RHO = 1 -> scale a little bit!" << std::endl;
				rho -= std::numeric_limits<double>::epsilon() * 10;
			}

			// exponential waiting time
			double dt = -log(rho) / R_tot;
			if (std::isnan(dt) || dt < 1.0e-15) {
				std::cout << "dt!->rho:  " << rho << " total rate: " << R_tot << "-log(rho)" << -log(rho) << std::endl;
			}

			// Strang first step
			//gp.timestep(dt / 2.0, dTdt, 0);
			t += (dt / 2.0); // First half-step
			if (t > end_time)
				break;

			gp.timestep_temp_grad(t, dt / 2.0, {0.0, 0.0});
			pp.volume_expansion(dt / 2.0, gp);

			// Strang second step
			g_si += pp.timestep(dt, gp, cnt);
			//gp.timestep(dt, 0, 0, { -g_si,0.0 });
			gp.timestep_temp_grad(t, dt / 2.0, {-g_si, 0.0});

			// Strang third step
			//gp.timestep(dt / 2.0, dTdt, 0);
			t += (dt / 2.0); // Second half-step
			if (t > end_time)
				break;

			gp.timestep_temp_grad(t, dt / 2.0, {0.0, 0.0});
			pp.volume_expansion(dt / 2.0, gp);

			//t += dt;
			iter++;

			if (counter_trigger(iter, PRINT_EVERY)) {

				clock.stop();

				double T = gp.get_T();
				double S = gp.get_S("Si");

				std::cout
					<< "time: "			<< t << '\t'
					<< "Temp: "				<< gp.get_T() << '\t'
					<< "SuperSat: "				<< gp.get_S("Si") << '\t'
					<< "Nucl_Rate: " << cnt.nucleation_rate(T, S) << '\t'
					<< "Gas-Phase: "		<< gp.get_n("Si") << '\t'
					<< "Stable_C_size: "		<< cnt.stable_cluster_size(T, S) << '\t'
					<< "AVG_Part_per_Agg: "	<< pp.get_mean_particles_number() << '\t'
					<< "AVG_Sint_level: "	<< pp.get_mean_sintering_level() << '\t'
					<< "AVG Sperical_Diam: " << pp.get_aggregates_mean_spherical_diameter() << '\t'
					<< "Agg_Number: "		<< pp.get_aggregates_number() << '\t'
					<< "Agg_density: "	<< pp.get_aggregates_density() << '\t'
					<< "Volume: "			<< pp.get_volume() << '\t'
					<< "Fractal_Dim: "			<< pp.get_mean_fractal_dimension() << '\t'
					<< "ts_exec_time: "	<< clock.interval() / PRINT_EVERY << std::endl;

				// Print file

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

		plot_data.close();

		// Anyway at the end save aggregates data for PSD
		std::valarray<double> particles_sizes = pp.get_particles_sizes();
		std::valarray<double> aggregates_sizes = pp.get_aggregates_sizes();

		// print particles sizes
		for (int i = 0; i < particles_sizes.size(); i++) {
			part_sizes_file << particles_sizes[i] << " ";
			global_psd_primary << particles_sizes[i] << " ";
		}
		part_sizes_file << std::endl;
		// print aggregates sizes
		for (int i = 0; i < aggregates_sizes.size(); i++) {
			agg_sizes_file << aggregates_sizes[i] << " ";
			global_psd_aggregates << aggregates_sizes[i] << " ";
		}
		agg_sizes_file << std::endl;

		part_sizes_file.close();
		agg_sizes_file.close();


	} // End of streamlines cycle

	// close global PSD files
	global_psd_aggregates.close();
	global_psd_primary.close();

	system("PAUSE");

	return;

}
