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

#include "linked_langevin_simulation.h"

LinkedLangevinSimulation::LinkedLangevinSimulation(raw_configuration_data _xml_data) :
	LinkedSimulation(_xml_data),
	V(_xml_data.volume),
	max_particles(_xml_data.max_part),
	min_particles(_xml_data.min_part),
	VTK_PATH(_xml_data.vtk_save),
	SAVE_VTK(_xml_data.SAVE_VTK_STEPS),
	max_dt(_xml_data.max_time_step){

	// Print Simulation data
	std::cout << "Volume [m3]: " << V << std::endl;
	std::cout << "Max Particles [#]: " << max_particles << std::endl;
	std::cout << "Min Particles [#]: " << min_particles << std::endl;
	std::cout << "VTK Files Saving Path: " << VTK_PATH << std::endl;
	std::cout << "VTK Files Saving Frequency " << SAVE_VTK << std::endl;
	std::cout << "max dt allowed [sec]: " << max_dt << std::endl;

}

void LinkedLangevinSimulation::run_simulation() {
	WallClock clock;

	// Declare the path for the XML file
	//std::string xml_path = "streams/NanoDomeStreamlines.xml";
	std::string xml_path = this->STREAMS_PATH;

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Global PSD from all streamlines (primary particles)
	std::ofstream global_psd_primary;
	std::string global_psd_primary_file = this->SAVE_PATH + "global_primary_PSD.dat";
	global_psd_primary.open(global_psd_primary_file);

	// Global PSD from all streamlines (aggregates)
	std::ofstream global_psd_aggregates;
	std::string global_psd_aggregates_file = this->SAVE_PATH + "global_aggregates_PSD.dat";
	global_psd_aggregates.open(global_psd_aggregates_file);

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	/// Execute a simulation for each streamline in the read file

	for (int s = 0; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;

		// Creates Streamline folder
		std::string stream_folder = this->SAVE_VTK + "Streamline" + std::to_string(s) + "/";

		std::experimental::filesystem::create_directory(stream_folder);

		// Create the species list
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
		double dt = 1.0e-8;
		// Temperature
		double T = 0.0;

		GasPhaseLink gp(l_streams_xml[s].get_splineT(),
			l_streams_xml[s].get_splineP(),
			gas_species, l_sc,
			l_streams_xml[s].get_minTime());


		// setup the nucleation theory we want to use
		ClassicalNucleationTheory cnt(gas_species[0]);

		// Starting pressure 
		double V_start = this->V;

		// setup the particle phase
		ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> pp(V_start);

		double elapsed_time = l_streams_xml[s].get_minTime(); // Simulation time
		double end_time = l_streams_xml[s].get_maxTime(); // Simulation end Time

		int iterations = 0; // Iterations
		double snap_count = 0.0; // Counter for saving VTK file

		const int PRINT_EVERY = this->PRINT_STEPS;
		double SAVE_SNAPSHOT = this->SAVE_VTK;
		const int SAVE_EVERY = this->SAVE_STEPS;

		std::string dir = this->SAVE_PATH;
		std::ofstream plot_data;
		std::string filename = dir + "Langevin_plot" + std::to_string(s) + "Streamline.dat";
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

		// PSD DATA
		std::ofstream part_sizes_file;
		std::string part_size_filename = dir + "particles_sizes"+ std::to_string(s) +"Streamline.dat";
		part_sizes_file.open(part_size_filename);
		// Aggregates diameters
		std::ofstream agg_sizes_file;
		std::string agg_size_filename = dir + "aggregates_sizes" + std::to_string(s) + "Streamline.dat";
		agg_sizes_file.open(agg_size_filename);

		// Simulation Main Cycle
		while (elapsed_time < end_time) {

			// check the smallest particle in the system
			double d_min = pp.get_particles_smallest_diameter();

			// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
			double dt_max_lang = d_min*gas_species[0].get_bulk_density(gp.get_T()) / gp.get_gas_flux();
			// calculate dt max to have v*dt < d/2 for the smallest particle
			double dt_max_coll = sqrt(M_PI*gas_species[0].get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

			if (pp.get_aggregates_number() >= 1) {
				dt = std::min(dt_max_coll, dt_max_lang);
				if (dt > this->max_dt) dt = 1.0e-8;
			}

			T = gp.get_T();

			// particle phase timestep
			double g_si = pp.timestep_lc(dt, gp, cnt, T);

			// gas phase time step
			gp.timestep(elapsed_time, dt);

			// Update elapsed time and iterations
			elapsed_time += dt;
			snap_count += dt;
			iterations++;

			if (counter_trigger(iterations, PRINT_EVERY)) {

				clock.stop();

				double T = gp.get_T();
				double S = gp.get_S("Si");

				std::cout
					<< "time: "			<< elapsed_time << '\t'
					<< "Temp: "				<< gp.get_T() << '\t'
					<< "SuperSat: "				<< gp.get_S("Si") << '\t'
					<< "Nucl_Rate: " << cnt.nucleation_rate(T, S) << '\t'
					<< "Gas-Phase: "		<< gp.get_n("Si") << '\t'
					<< "Stable_C_size: "		<< cnt.stable_cluster_size(T, S) << '\t'
					<< "AVG_Part_per_Agg: "	<< pp.get_mean_particles_number() << '\t'
					<< "AVG_Sint_level: "	<< pp.get_mean_sintering_level() << '\t'
					<< "AVG Sperical_Diam: " << pp.get_aggregates_mean_spherical_diameter() << '\t'
					<< "Agg_Numbe: "		<< pp.get_aggregates_number() << '\t'
					<< "Agg_density: "	<< pp.get_aggregates_density() << '\t'
					<< "Volume: "			<< pp.get_volume() << '\t'
					<< "ts: "				<< dt << '\t'
					<< "Fractal_Dim: "			<< pp.get_mean_fractal_dimension() << '\t'
					<< "ts_exec_time: "	<<clock.interval() / PRINT_EVERY << std::endl;

				// Print file

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

				// Print file

				clock.start();
			}

			//// Print snapshot
			//if (iterations%PRINT_STEP == 0 && pp.get_aggregates_number() > 0) {
			//	std::cout << elapsed_time << '\t'
			//		<< dt_max_coll << '\t'
			//		<< dt_max_lang << '\t'
			//		<< "T[K]: " << gp.get_T() << '\t'
			//		<< "V[m^3]" << pp.get_volume() << '\t'
			//		//<< J_si << '\t'
			//		//<< ns << '\t'
			//		//<< j_si << '\t'
			//		<< "|N|: " << pp.get_aggregates_number()
			//		//<< "|N|: " << pp.get_aggregates_cardinality()
			//		<< std::endl;
			//}

			// Save VTK
			if (snap_count >= SAVE_SNAPSHOT && pp.get_aggregates_number() > 0) {
				pp.save_vtk(iterations, stream_folder);
				snap_count = 0.0;
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

	}

	// close global PSD files
	global_psd_aggregates.close();
	global_psd_primary.close();

	return;

}
