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

#include "utilities.h"
#include "clock.h"
#include "ndm_random.h"

#include "gasphase/gasphasecv.h"
#include "gasphase/gasphaselink.h"


#include "fmmajorantkernel.h"

#include "splinelinear.h"

// OpenMP
#include <omp.h>


// I/O
#include "i_o.h"
#include "streamlinefileXML.h"
#include "streamlinefileJSON.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <memory>
#include <iostream>
#include <fstream>

/// AUXILIARY FUNCTIONS
double Temp_gradient(double _t, double _T) {

	double cycle = 1.6;
	double period = 0.4;
	int n_cycle = std::floor(_t / cycle);
	double sub_cycle = _t - n_cycle*cycle;
	int index = std::floor(sub_cycle / period);

	if (index == 0)
		return -1.0e4;
	else if (index == 1)
		return 0.0;
	else if (index == 2)
		return 1.0e4;
	else if (index == 3)
		return 0.0;

	/*
	double cycle = 0.0016;
	double period = 0.0004;
	int n_cycle = std::floor(_t / cycle);
	double sub_cycle = _t - n_cycle*cycle;
	int index = std::floor(sub_cycle / period);

	if (index == 0)
		return -1.0e7;
	else if (index == 1)
		return 1.0e7;
	else if (index == 2)
		return -1.0e7;
	else if (index == 3)
		return 1.0e7;
		*/

	/*
	double cycle = 0.016;
	double period = 0.004;
	int n_cycle = std::floor(_t / cycle);
	double sub_cycle = _t - n_cycle*cycle;
	int index = std::floor(sub_cycle / period);

	if (index == 0)
		return -1.0e6;
	else if (index == 1)
		return 0.0;
	else if (index == 2)
		return 1.0e6;
	else if (index == 3)
		return 0.0;
		*/

	/*double cycle = 0.00108;
	int n_cycle = _t / cycle;
	int sub_cycle = _t % cycle;


	if (_t / step*c < 1.0)
		return -1.0e7;
	else if (_t / step*c >= 1.0 && _t / step*c < 2.0)
		return 0.0;
	else if (_t / step*c >= 2.0 && _t / step*c < 3.0)
		return 1.0e7;
	else if (_t / step*c >= 3.0 && _t / step*c < 4.0)
		return 0.0;*/
	
	/*double complete = 0.00108;

	if (_t <= 0.00027 && _t  > 0.0)
		return -1.0e7;
	else if (_t  >= 0.00027 && _t  < 0.00054)
		return 0.0;
	else if (_t  >= 0.00054 && _t < 0.00081)
		return 1.0e7;
	else if (_t  >= 0.00081 && _t < 0.00108)
		return 0.0;
	else if (_t>= complete && _t <=  0.00027 + complete)
		return -1.0e7;
	else if (_t >= 0.00027 + complete && _t  < 0.00054 + complete)
		return 0.0;
	else if (_t >= 0.00054 + complete && _t < 0.00081 + complete)
		return 1.0e7;
	else if (_t >= 0.00081 + complete && _t < 0.00108 + complete)
		return 0.0;*/

}

double Temp_gradient(double _t, double _T, double _t_ref) {

	double cycle = 1.6;
	double period = 0.4;
	int n_cycle = std::floor(_t / cycle);
	double sub_cycle = _t - n_cycle*cycle;
	int index = std::floor(sub_cycle / period);

	if (index == 0)
		return -1.0e4;
	else if (index == 1)
		return 0.0;
	else if (index == 2 && _T < _t_ref)
		return 1.0e4;
	else if(index == 2 && _T >= _t_ref)
		return 0.0;
	else if (index == 3)
		return 0.0;

	/*
	double cycle = 0.0016;
	double period = 0.0004;
	int n_cycle = std::floor(_t / cycle);
	double sub_cycle = _t - n_cycle*cycle;
	int index = std::floor(sub_cycle / period);

	if (index == 0)
	return -1.0e7;
	else if (index == 1)
	return 1.0e7;
	else if (index == 2)
	return -1.0e7;
	else if (index == 3)
	return 1.0e7;
	*/

	/*
	double cycle = 0.016;
	double period = 0.004;
	int n_cycle = std::floor(_t / cycle);
	double sub_cycle = _t - n_cycle*cycle;
	int index = std::floor(sub_cycle / period);

	if (index == 0)
	return -1.0e6;
	else if (index == 1)
	return 0.0;
	else if (index == 2)
	return 1.0e6;
	else if (index == 3)
	return 0.0;
	*/

	/*double cycle = 0.00108;
	int n_cycle = _t / cycle;
	int sub_cycle = _t % cycle;


	if (_t / step*c < 1.0)
	return -1.0e7;
	else if (_t / step*c >= 1.0 && _t / step*c < 2.0)
	return 0.0;
	else if (_t / step*c >= 2.0 && _t / step*c < 3.0)
	return 1.0e7;
	else if (_t / step*c >= 3.0 && _t / step*c < 4.0)
	return 0.0;*/

	/*double complete = 0.00108;

	if (_t <= 0.00027 && _t  > 0.0)
	return -1.0e7;
	else if (_t  >= 0.00027 && _t  < 0.00054)
	return 0.0;
	else if (_t  >= 0.00054 && _t < 0.00081)
	return 1.0e7;
	else if (_t  >= 0.00081 && _t < 0.00108)
	return 0.0;
	else if (_t>= complete && _t <=  0.00027 + complete)
	return -1.0e7;
	else if (_t >= 0.00027 + complete && _t  < 0.00054 + complete)
	return 0.0;
	else if (_t >= 0.00054 + complete && _t < 0.00081 + complete)
	return 1.0e7;
	else if (_t >= 0.00081 + complete && _t < 0.00108 + complete)
	return 0.0;*/

}



// Moments method
#include "moments/momentmodelfriedlander.h"
#include "moments/momentmodelpratsinis.h"

void moment() {

    Species si("Si");
    Species ar("Ar");

    // create the vector of species for the gas phase
    std::vector<Species> gas_species = {si,ar};

    // simulation timestep
    const double dt = 1.0e-8;
	std::string s_dt = "";

    // gas phase info
    double p = 1.01e5; // pressure [Pa]
    //double dTdt = 0.0; // temperature gradient [K/s]
	double dTdt = 1.0e7;
	//std::string s_dTdt = "0.0";

    double T_start = 2000.0;

    // set pressure, temperature, species and initial molar concentration
    //GasPhaseCV gp(p, T_start, gas_species, {0.25, 0.75});
	GasPhaseCV gp(p, T_start, gas_species, { 0.0, 1.0 });

    // setup the nucleation theory we want to use
    ClassicalNucleationTheory cnt(si);

    // set the moment method with the condensing species
    //MomentModelFriedlander mm(si);
    //MomentModelPratsinis mm(si);
	MomentModelPratsinis mm(si, 5.17e19, 3.82e-6, 3.58e-31);

    double t = 0.;
    int iter = 0;
	
	const int PRINT_EVERY = 1000;

	const int SAVE_EVERY = 1;

	// Clean lognormal
	std::string lognormal_path = "Lognormal" + s_dt + ".dat";
	std::ofstream lognormal_f;
	lognormal_f.open(lognormal_path);
	lognormal_f.clear();
	lognormal_f.close();

	std::ofstream plot_data;
	std::string filename = "MOMENTS_plot" + s_dt + ".dat";
	plot_data.open(filename);
	// Print headlines
	/*plot_data << "Time[sec]" << '\t'
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
		<< std::endl;*/

	

    // loop over timesteps until the final temperature
    while(t < 1.5) {

        double T = gp.get_T();
        double ns = gp.get_n("Si");
        double S = ns/si.n_sat(T);
        double J = cnt.nucleation_rate(T,S);
        double j = cnt.stable_cluster_size(T,S);

		if (T >= 3500) break;

		//dTdt = Temp_gradient(t, T, 1000);
		//dTdt = Temp_gradient(t, 2);

        // moment method timestep and species consumption retrieval
        double g_si = mm.timestep(dt,gp,cnt);

        // updating the gas phase
        gp.timestep(dt, dTdt, 0, {-g_si, 0.0});

        ++iter;
        t+=dt;

		if (counter_trigger(iter, PRINT_EVERY)) {
			//gp.print();

			std::cout
				<< "Time= "<< t << '\t'
				<< "Temp= " << T << '\t'
				<< "Sat= " << S << '\t'
				<< "Nucl Rate= " << J << '\t'
				<< "GP[mon]= "<< ns << '\t'
				<< "s.c.size[mon]= "<< j << '\t'
				<< "Density= " << mm.get_density() << '\t'
				<< "Mean Diam= " << mm.get_mean_diameter() << '\t'
				<< "GP_cond= "<<mm.get_cond_term() << '\t'
				<< "M1= "<<mm.get_M1()<<'\t'
				<< "M2= " << mm.get_M2() << '\t'
				<< std::endl;
		}

		// File savings
		if (counter_trigger(iter, SAVE_EVERY)) {

			// save lognormal values
			mm.get_lognormal_val();
			// save simulation data
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
				<< 0.0 <<  '\t' // dummy interval time
				<< mm.get_cond_term() << '\t'
				<< mm.get_M1() << '\t'
				<< mm.get_M2()
				<< std::endl;
		}
    }

	plot_data.close();

	system("PAUSE");
}

void moment_link() {

    WallClock clock;

    // for testing purposes
    // Declare the path for the XML file
    //std::string xml_path = "streams/NanoDomeStreamlines_brux3.xml";
	std::string xml_path = "streams/NanoDomeStreamlines_test_2018_01_12.xml";
	//std::string xml_path = "streams/NanoDomeStreamlinesMarcus2.xml"; 
    // Declare the path for the JSON file
    //std::string json_path = "streams/NanoDomeStreamlines.json";

    /*---------------------- XML TESTING -----------------------------*/

    // Create the vector containing all the streamlines in the file
    std::vector<streamline> l_streams_xml;

    // Inizialize the data structure to manage the XML file
    streamlinefileXML s_xml(xml_path, l_streams_xml);

    // Create the streamline
    s_xml.read_Streamlines();

    Species si("Si");
    Species ar("Ar");

    // create the vector of species for the gas phase
    std::vector<Species> gas_species = { si,ar };

    // Main cycle: a simulation for each streamline extracted from the CFD

    for(int s = 0; s < l_streams_xml.size(); s++){

        std::cout << "STREAMLINE: " << s << std::endl;

        // Create splines interpolating the samples extracted
        l_streams_xml[s].create_splines();

        std::vector<Spline*> l_sc;
        l_sc.push_back(l_streams_xml[s].get_splineMF());

        // Dummy Spline for Ar
        SplineLinear *sC_ar = new SplineLinear();
        l_sc.push_back(sC_ar);

        // simulation timestep
        double dt = 1e-9;

        GasPhaseLink gp(l_streams_xml[s].get_splineT(),
                        l_streams_xml[s].get_splineP(),
                        gas_species, l_sc,
                        l_streams_xml[s].get_minTime());


        // setup the nucleation theory we want to use
        ClassicalNucleationTheory cnt(si);

        // set the moment method with the condensing species
        //MomentModelFriedlander mm(si);
        MomentModelPratsinis mm(si);

        double t = l_streams_xml[s].get_minTime();
		//double t = 0.1;

        // simulation end time
        double end_time = l_streams_xml[s].get_maxTime();

        int iter = 0;

        int PRINT_EVERY = 10000;

        const int SAVE_EVERY = 100;
        std::ofstream plot_data;
		std::ofstream plotM1;
        std::string filename = "MOMENTS_plot_" + std::to_string(s) + "_Streamline.dat";
        plot_data.open(filename);

		// Clean lognormal
		std::string lognormal_path = "Lognormal_" + std::to_string(s) + ".dat";
		std::ofstream lognormal_f;
		lognormal_f.open(lognormal_path);
		lognormal_f.clear();
		lognormal_f.close();

		// Print headlines
		/*plot_data 
			<< "1)Time[sec]" << '\t'
			<< "2)Temp[K]" << '\t'
			<< "3)SuperSaturation ratio" << '\t'
			<< "4)Nucl Rate" << '\t'
			<< "5)Species # density" << '\t'
			<< "6)Stable cluster size[m]" << '\t'
			<< "7)AVG Part Num[#]: " << '\t'
			<< "8)Sint level[%]" << '\t'
			<< "9)AVG diameter[m]" << '\t'
			<< "10)Agg. #[#]" << '\t'
			<< "11)Agg density[#/m3]" << '\t'
			<< "12)Volume[m3]" << '\t'
			<< "13)AVG fract dim" << '\t'
			<< "14)ts exec time"
			<< std::endl;*/

        // loop over timesteps until the final temperature
        while (t < end_time) {

            double T = gp.get_T();
            double ns = gp.get_n("Si");
            double S = ns / si.n_sat(T);
            double J = cnt.nucleation_rate(T, S);
            double j = cnt.stable_cluster_size(T, S);

            // moment method timestep and species consumption retrieval
            double g_si = mm.timestep(dt, gp, cnt, iter);

            // updating the gas phase
            //gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });
            gp.timestep(t, dt);


			/*if (t > 1.6665) {
				dt = 1.0e-10;
				
			}

			if (iter > 18272595)
				PRINT_EVERY = 1;
			*/

            ++iter;
            t += dt;

            // Step Print
            if (counter_trigger(iter, PRINT_EVERY)) {
                //gp.print();
				int std_prec = std::cout.precision();
				std::cout.precision(15);
				std::cout
					<< t << '\t';
				std::cout.precision(std_prec);
				std::cout
					<< iter << '\t'
					<< T << '\t'
					<< S << '\t'
					<< J << '\t'
					<< ns << '\t'
					<< j << '\t'
					<< mm.get_density() << '\t'
					<< mm.get_mean_diameter() << '\t'
					<< gp.get_molar_f()
					<< std::endl;
            }

            // File savings
            if (counter_trigger(iter, SAVE_EVERY)) {

				mm.get_lognormal_val(s);

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

		mm.get_lognormal_val(s);

		double T = gp.get_T();
		double ns = gp.get_n("Si");
		double S = ns / si.n_sat(T);
		double J = cnt.nucleation_rate(T, S);
		double j = cnt.stable_cluster_size(T, S);

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
		plot_data.close();
    }

    system("PAUSE");

}

void moment_link_temp() {

	WallClock clock;

	// for testing purposes
	// Declare the path for the XML file
	//std::string xml_path = "streams/NanoDomeStreamlines.xml";
	//std::string xml_path = "streams/NanoDomeStreamlines_Marcus_test.xml";
	//std::string xml_path = "streams/NanoDomeStreamlinesMarcus2.xml"; 
	std::string xml_path = "streams/NanoDomeStreamlines_test_2018_01_12.xml";
	// Declare the path for the JSON file
	//std::string json_path = "streams/NanoDomeStreamlines.json";

	/*---------------------- XML TESTING -----------------------------*/

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si,ar };

	// Main cycle: a simulation for each streamline extracted from the CFD

	for (int s = 3; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;

		// Create splines interpolating the samples extracted
		l_streams_xml[s].create_splines();
		
		// simulation timestep
		const double dt = 1e-8;
		// streamline start time
		//double start_time = l_streams_xml[s].get_minTime();
		double start_time = 0.12;
		// starting molar fraction value
		double start_mf = l_streams_xml[s].get_splineMF()->get_v(start_time);

		GasPhaseLink gp(l_streams_xml[s].get_splineT(),
			l_streams_xml[s].get_splineP(),
			gas_species, start_time, start_mf);


		// setup the nucleation theory we want to use
		ClassicalNucleationTheory cnt(si);

		// set the moment method with the condensing species
		//MomentModelFriedlander mm(si);
		MomentModelPratsinis mm(si);

		double t = start_time;

		// simulation end time
		double end_time = l_streams_xml[s].get_maxTime();

		int iter = 0;

		const int PRINT_EVERY = 10000;

		const int SAVE_EVERY = 100;
		std::ofstream plot_data;
		std::ofstream plotM1;
		std::string filename = "MOMENTS_plot" + std::to_string(s) + "Streamline.dat";
		plot_data.open(filename);
		// Print headlines
		/*plot_data << "Time[sec]" << '\t'
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
			<< std::endl;*/

		// loop over timesteps until the final temperature
		while (t < end_time) {

			double T = gp.get_T();
			double ns = gp.get_n("Si");
			double S = ns / si.n_sat(T);
			double J = cnt.nucleation_rate(T, S);
			double j = cnt.stable_cluster_size(T, S);

			// moment method timestep and species consumption retrieval
			double g_si = mm.timestep(dt, gp, cnt);

			// updating the gas phase
			gp.timestep_temp_grad(t, dt, {-g_si, 0.0});

			++iter;
			t += dt;

			// Step Print
			if (counter_trigger(iter, PRINT_EVERY)) {
				//gp.print();

				std::cout
					<< t << '\t'
					<< T << '\t'
					<< S << '\t'
					<< J << '\t'
					<< ns << '\t'
					<< j << '\t'
					<< mm.get_density() << '\t'
					<< mm.get_mean_diameter() << std::endl;
			}

			// File savings
			if (counter_trigger(iter, SAVE_EVERY)) {

				mm.get_lognormal_val(s);

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

		mm.get_lognormal_val(s);
		plot_data.close();
	}

	system("PAUSE");
}

/*---------- TEST PARALLEL IMPLEMENTATION -----------------------------*/
std::vector<std::vector<double>> moment_link_thd(streamline& s, int s_idx) {

    WallClock clock;

    Species si("Si");
    Species ar("Ar");

    // create the vector of species for the gas phase
    std::vector<Species> gas_species = { si,ar };

    // Create Spline for Temperature
    SplineLinear *sT = new SplineLinear(s.get_Time(), s.get_Temp());

    // Create Spline for Pressure
    SplineLinear *sP = new SplineLinear(s.get_Time(), s.get_Press());

    // Create list of Spline for species concentration
    SplineLinear *sC_si = new SplineLinear(s.get_Time(), s.get_Molar());
    //std::vector<Spline*> l_sc; l_sc.push_back(sC_si);

    // Create list of Spline for species molar fraction
    SplineLinear *sC_f_si = new SplineLinear(s.get_Time(), s.get_Molar_Fraction());
    std::vector<Spline*> l_sc; l_sc.push_back(sC_f_si);

    // Dummy Spline for Ar
    SplineLinear *sC_ar = new SplineLinear();
    l_sc.push_back(sC_ar);

    // simulation timestep
    const double dt = 1e-8;

    // gas phase info
    //double p = 1.01e5; // pressure [Pa]
    //double dTdt = -1e7; // temperature gradient [K/s]

    //double T_start = 3000.0;
    double T_start = sT->get_v(sT->get_x_min());

    // set pressure, temperature, species and initial molar concentration
    //GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75 });
    GasPhaseLink gp(sT, sP, gas_species, l_sc, sT->get_x_min());

    // setup the nucleation theory we want to use
    ClassicalNucleationTheory cnt(si);

    // set the moment method with the condensing species
    //MomentModelFriedlander mm(si);
    MomentModelPratsinis mm(si);

    double t = sT->get_x_min();
    double end_time = sT->get_x_max();

    int iter = 0;

    const int PRINT_EVERY = 500;

    const int SAVE_EVERY = 30;
    std::vector<std::vector<double>> plot_data;

    // loop over timesteps until the final temperature
    while (t < end_time) {

        double T = gp.get_T();
        double ns = gp.get_n("Si");
        double S = ns / si.n_sat(T);
        double J = cnt.nucleation_rate(T, S);
        double j = cnt.stable_cluster_size(T, S);

        // moment method timestep and species consumption retrieval
        double g_si = mm.timestep(dt, gp, cnt);

        // updating the gas phase
        //gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

        gp.timestep(t, dt);

        ++iter;
        t += dt;

        // Terminal Print
        //if (counter_trigger(iter, PRINT_EVERY)) {
        //	//gp.print();

        //	std::cout << t << '\t' << T << '\t' << S << '\t' << J << '\t' << ns << '\t'
        //		<< j << '\t' << mm.get_density() << '\t' << mm.get_mean_diameter() << std::endl;
        //}

        // plot entry;
        if (counter_trigger(iter, SAVE_EVERY)) {
            std::vector<double> buffer;
            buffer.push_back(t); buffer.push_back(T); buffer.push_back(S); buffer.push_back(J);
            buffer.push_back(ns); buffer.push_back(j); buffer.push_back(mm.get_density());
            buffer.push_back(mm.get_mean_diameter());
            /* << t << '\t' << T << '\t' << S << '\t' << J << '\t' << ns << '\t'
                << j << '\t' << mm.get_density() << '\t' << mm.get_mean_diameter() << std::endl;*/
            plot_data.push_back(buffer);
        }

    }

    return plot_data;
}

void moment_parallel() {

    // READ STREAMLINES FROM THE XML FILE

    // for testing purposes
    // Declare the path for the XML file
    std::string xml_path = "streams/NanoDomeStreamlines.xml";
    // Declare the path for the JSON file
    //std::string json_path = "streams/NanoDomeStreamlines.json";

    /*---------------------- XML TESTING -----------------------------*/

    // Create the vector containing all the streamlines in the file
    std::vector<streamline> l_streams_xml;

    // Inizialize the data structure to manage the XML file
    streamlinefileXML s_xml(xml_path, l_streams_xml);

    // Create the streamline
    s_xml.read_Streamlines();

    // vector for simulations data
    std::vector<std::vector<std::vector<double>>> results;
    for (int i = 0; i < l_streams_xml.size(); i++) {
        std::vector<std::vector<double>> blank;
        results.push_back(blank);
    }

    // Parallel execution of the simulations
    WallClock clock;
    omp_set_num_threads(8);
    int s;
    std::vector<std::vector<double>> p_res;

    clock.start();

#pragma omp parallel for private(s, p_res)
    for (s = 0; s < l_streams_xml.size(); s++) {

        p_res = moment_link_thd(l_streams_xml[s], s);
        results[s] = p_res;

    }

    clock.stop();
    std::cout<<"ELAPSED TIME: "<< clock.interval()<<std::endl;

    // Print Plot Files
    for (int i = 0; i < results.size(); i++) {

        std::ofstream file;
        std::string filename = "MOMENTS_plot" + std::to_string(i) + "Streamline.dat";
        file.open(filename);
        for (int t = 0; t < results[i].size(); t++) {
            for (int v = 0; v < results[i][v].size(); v++) {
                file << results[i][t][v] << "\t";
            }
            file << std::endl;
        }

        file.close();
    }

    system("PAUSE");

}
/*---------- TEST PARALLEL IMPLEMENTATION -----------------------------*/

// PBM model
#include "particlephase/pbmfractalparticlephase.h"
#include "aggregate/pbmaggregate.h"

void particle() {

	/*	On a multivariate population balance model to describe the structure and composition of silica nanoparticles.
		Shraddha Shekara, William J. Menza, Alastair J. Smith, Markus Kraft, Wolfgang Wagnerb
		Computers and Chemical Engineering	*/

    WallClock clock;

    Species si("Si");
    Species ar("Ar");

    // create the vector of species for the gas phase
    std::vector<Species> gas_species = {si, ar};

    // gas phase info
    double p = 1.01e5; // pressure [Pa]
    double dTdt = -1.0e7; // temperature gradient [K/s]

    double T_start = 3000.0;

    // set pressure, temperature, species and initial molar concentration
    GasPhaseCV gp(p, T_start, gas_species, {0.01, 0.99});

    // setup the nucleation theory we want to use
    ClassicalNucleationTheory cnt(si);

    // setup the particle phase we want to use
    PBMFractalParticlePhase<PBMAggregate<Particle>> pp(1.61, 5.0e-16);

	// set max and min number of aggregates
	pp.set_max_aggregates(2000);
	pp.set_min_aggregates(1990);

    double t = 0.0;
    int iter = 0;

    const int PRINT_EVERY = 1000;
	const int PRINT_HEADLINE = PRINT_EVERY * 5;


	const int SAVE_EVERY = 1000;
	std::ofstream plot_data;
	std::string filename = "PBM_plot.dat";
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

	// PSD Data files
	const double PSD_DATA = 1.0e-5;
	double save_step = 0.0;
	// Prinary particles
	std::ofstream part_sizes_file;
	std::string part_size_filename = "particles_sizes.dat";
	part_sizes_file.open(part_size_filename);
	// Aggregates diameters
	std::ofstream agg_sizes_file;
	std::string agg_size_filename = "aggregates_sizes.dat";
	agg_sizes_file.open(agg_size_filename);

    clock.start();

    // loop over timesteps
    while(t < 0.002) {

		if (gp.get_T() < 300) 
		{ dTdt = 0; }

        // species source term for the gas phase
        double g_si = 0.0;

        // calculate the timestep using an exponential waiting time
        double R_tot = pp.get_total_processes_rate(gp,cnt);
        double rho = ndm::uniform_double_distr(ndm::rand_gen);

        // exponential waiting time
        double dt = -log(rho)/R_tot;

        // Strang first step
        gp.timestep(dt/2.0,dTdt,0);
        pp.volume_expansion(dt/2.0,gp);

        // Strang second step
        g_si += pp.timestep(dt,gp,cnt);
        gp.timestep(dt,0,0,{-g_si,0.0});

        // Strang third step
        gp.timestep(dt/2.0,dTdt,0);
        pp.volume_expansion(dt/2.0,gp);

        t += dt;
		save_step += dt;
        iter++;

        if(counter_trigger(iter,PRINT_EVERY)) {

            clock.stop();

            double T = gp.get_T();
            double S = gp.get_S("Si");

            std::cout << t << '\t'									// time
                      << gp.get_T() << '\t'							// temperature
                      << gp.get_S("Si") << '\t'						// supersaturation (S)
                      << cnt.nucleation_rate(T, S) << '\t'			// J
                      << gp.get_n("Si") << '\t'						// ns
                      << cnt.stable_cluster_size(T, S) << '\t'		// j
                      << pp.get_mean_particles_number() << '\t'		// N_m
                      << pp.get_mean_sintering_level() << '\t'		//
                      << pp.get_aggregates_mean_spherical_diameter() << '\t' 
                      << pp.get_aggregates_number() << '\t'
                      << pp.get_aggregates_density() << '\t'
					  << pp.get_volume() << '\t'
					  << pp.get_mean_fractal_dimension() << '\t'
                      << clock.interval()/PRINT_EVERY << std::endl;

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
			
			clock.start();
		}

		// Save particles sizes for PSD
		//if (save_step>=PSD_DATA) {

		//	clock.stop();

		//	std::valarray<double> particles_sizes = pp.get_particles_sizes();
		//	std::valarray<double> aggregates_sizes = pp.get_aggregates_sizes();

		//	// print particles sizes
		//	for (int i = 0; i < particles_sizes.size(); i++) {
		//		part_sizes_file << particles_sizes[i] << " ";
		//	}
		//	part_sizes_file << std::endl;
		//	// print aggregates sizes
		//	for (int i = 0; i < aggregates_sizes.size(); i++) {
		//		agg_sizes_file << aggregates_sizes[i] << " ";
		//	}
		//	agg_sizes_file << std::endl;

		//	save_step = 0.0;

		//	clock.start();
		//}
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

}

void particle_link() {

    WallClock clock;

    // for testing purposes
    // Declare the path for the XML file
    //std::string xml_path = "streams/NanoDomeStreamlines_brux3.xml";
	std::string xml_path = "streams/NanoDomeStreamlines_test_2018_01_12.xml";
	//std::string xml_path = "streams/TestStream.xml";
    // Declare the path for the JSON file
    //std::string json_path = "streams/NanoDomeStreamlines.json";

	// Global PSD from all streamlines (primary particles)
	std::ofstream global_psd_primary;
	std::string global_psd_primary_file = "global_primary_PSD.dat";
	global_psd_primary.open(global_psd_primary_file);

	// Global PSD from all streamlines (aggregates)
	std::ofstream global_psd_aggregates;
	std::string global_psd_aggregates_file = "global_aggregates_PSD.dat";
	global_psd_aggregates.open(global_psd_aggregates_file);

	// Create the vector containing all the streamlines in the file
    std::vector<streamline> l_streams_xml;

    // Inizialize the data structure to manage the XML file
    streamlinefileXML s_xml(xml_path, l_streams_xml);

    // Create the streamline
    s_xml.read_Streamlines();

    // FOR NOW NO SPECIES MANAGEMENT FROM THE STREAMLINE
    Species si("Si");
    Species ar("Ar");

    // create the vector of species for the gas phase
    std::vector<Species> gas_species = { si, ar };


    /// Execute a simulation for each streamline in the XML file

	// Streamlines loop
    for (int s = 0; s < l_streams_xml.size(); s++) {

        std::cout << "STREAMLINE: " << s << std::endl;

        // create splines with samples
        l_streams_xml[s].create_splines();

        // get spline for the Si molar fraction
        std::vector<Spline*> l_sc;
        l_sc.push_back(l_streams_xml[s].get_splineMF());

        // Dummy Spline for Ar molar fraction
        SplineLinear *sC_ar = new SplineLinear();
        l_sc.push_back(sC_ar);

        // set pressure, temperature, species and initial molar concentration
        GasPhaseLink gp(l_streams_xml[s].get_splineT(),
                        l_streams_xml[s].get_splineP(),
                        gas_species, l_sc, l_streams_xml[s].get_minTime());

        // setup the nucleation theory we want to use
        ClassicalNucleationTheory cnt(si);

        // setup the particle phase we want to use
        PBMFractalParticlePhase<PBMAggregate<Particle>> pp(1.6, 1.0e-9);

		// set max and min number of aggregates
		pp.set_max_aggregates(2000);
		pp.set_min_aggregates(1990);

        // get the streamline starting time
        double t = l_streams_xml[s].get_minTime();

        // get the end time for the streamline
        double end_time = l_streams_xml[s].get_maxTime();

        int iter = 0;

        const int PRINT_EVERY = 10000;
		const int PRINT_HEADLINE = PRINT_EVERY * 10;

		const int SAVE_EVERY = 2000;
		std::ofstream plot_data;
        std::string filename = "PBM_plot" + std::to_string(s) + "Streamline.dat";
        plot_data.open(filename);
		// Print headlines
		plot_data << "1)Time[sec]" << '\t'
			<< "2)Temp[K]" << '\t'
			<< "3)SuperSaturation ratio" << '\t'
			<< "4)Nucl Rate" << '\t'
			<< "5)Species # density" << '\t'
			<< "6)Stable cluster size[m]" << '\t'
			<< "7)AVG Part Num[#]: " << '\t'
			<< "8)Sint level[%]" << '\t'
			<< "9)AVG diameter[m]" << '\t'
			<< "10)Agg. #[#]" << '\t'
			<< "11)Agg density[#/m3]" << '\t'
			<< "12)Volume[m3]" << '\t'
			<< "13)AVG fract dim" << '\t'
			<< "14)ts exec time"
			<< std::endl;

		// PSD Data files
		const double PSD_DATA = 1.0e-5;
		double save_step = 0.0;
		// Prinary particles
		std::ofstream part_sizes_file;
		std::string part_size_filename = "particles_sizes" + std::to_string(s) + "Streamline.dat";
		part_sizes_file.open(part_size_filename);
		// Aggregates diameters
		std::ofstream agg_sizes_file;
		std::string agg_size_filename = "aggregates_sizes" + std::to_string(s) + "Streamline.dat";
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
			if(rho == 1.0) {
				std::cout << "RHO = 1 -> scale a little bit!" << std::endl;
				rho -= std::numeric_limits<double>::epsilon()*10;
			}
				

            // exponential waiting time
            double dt = -log(rho) / R_tot;
			// if not valid dt, print usefull data
			if (std::isnan(dt)) {

				std::ofstream DUMP;
				DUMP.open("DUMP.dat");

				int default_precision = std::cout.precision();
				std::cout.precision(15);

				std::cout << "dt: "<< dt <<"dt!->rho:  " << rho << " total rate: " << R_tot << "-log(rho)"<< -log(rho) <<std::endl;
				std::cout << "GAS PHASE:" << std::endl;
				std::cout << "T: " << gp.get_T() << '\t'
					<< "S: " << gp.get_S("Si") << '\t'
					<< "Volume: " << pp.get_volume() << 't'
					<< "Nucl. Rate: " << cnt.nucleation_rate(gp.get_T(), gp.get_S("Si")) << '\t'
					<< "Molar Concentration: " << gp.get_n("Si") << '\t'
					<< "Stable C size: " << cnt.stable_cluster_size(gp.get_T(), gp.get_S("Si")) << '\t' 
					<<" Condensation Rate: "<<cnt.condensation_rate(gp.get_T(), gp.get_S("Si")) << std::endl;
				std::cout << "Particle Phase: " << std::endl;
				std::cout << "Coagulation Rate: " << pp.coagulation_rate(gp.get_T())<<std::endl;
				std::cout << " K1: " << pp.K1 << "K2: " << pp.K2 << " K21: " << pp.K21 << " K22: " << pp.K22 << std::endl;

				// DUMP ALL ON FILE
				DUMP << "dt: " << dt << "dt!->rho:  " << rho << " total rate: " << R_tot << "-log(rho)" << -log(rho) << std::endl;
				DUMP << "GAS PHASE:" <<std::endl;
				DUMP << "T: " << gp.get_T() << '\t'
					<< "S: " << gp.get_S("Si") << '\t'
					<< "Volume: " << pp.get_volume() << 't'
					<< "Nucl. Rate: " << cnt.nucleation_rate(gp.get_T(), gp.get_S("Si")) << '\t'
					<< "Molar Concentration: " << gp.get_n("Si") << '\t'
					<< "Stable C size: " << cnt.stable_cluster_size(gp.get_T(), gp.get_S("Si")) << '\t'
					<< " Condensation Rate: " << cnt.condensation_rate(gp.get_T(), gp.get_S("Si")) << std::endl;
				DUMP << "Particle Phase: " << std::endl;
				DUMP << "Coagulation Rate: " << pp.coagulation_rate(gp.get_T()) << std::endl;
				DUMP << " K1: " << pp.K1 << "K2: " << pp.K2 << " K21: " << pp.K21 << " K22: " << pp.K22 << std::endl;

				// Dump all aggregates on file
				pp.dump_aggregates(DUMP);


				system("PAUSE");
				std::cout.precision(default_precision);
			}

            // Strang first step
            //gp.timestep(dt / 2.0, dTdt, 0);
            t += (dt / 2.0); // First half-step
            if (t > end_time)
                break;

            gp.timestep(t, dt / 2.0);
            pp.volume_expansion(dt / 2.0, gp);

            // Strang second step
            g_si += pp.timestep(dt, gp, cnt);
            //gp.timestep(dt, 0, 0, { -g_si,0.0 });
            gp.timestep(t, dt / 2.0);

            // Strang third step
            //gp.timestep(dt / 2.0, dTdt, 0);
            t += (dt / 2.0); // Second half-step
            if (t > end_time)
                break;

            gp.timestep(t, dt / 2.0);
            pp.volume_expansion(dt / 2.0, gp);

            //t += dt;
			save_step += dt;
            iter++;
			
			
			if (counter_trigger(iter, PRINT_EVERY)) {

                clock.stop();

                double T = gp.get_T();
                double S = gp.get_S("Si");

                std::cout 
					<< t << '\t'
                    << gp.get_T() << '\t'
                    << gp.get_S("Si") << '\t'
                    << cnt.nucleation_rate(T, S) << '\t'
                    << gp.get_n("Si") << '\t'
                    << cnt.stable_cluster_size(T, S) << '\t'
                    <<"AVG Part Num: "<< pp.get_mean_particles_number() << '\t'
                    <<"S_level: "<< pp.get_mean_sintering_level() << '\t'
                    << pp.get_aggregates_mean_spherical_diameter() << '\t'
                    << pp.get_aggregates_number() << '\t'
                    << pp.get_aggregates_density() << '\t'
                    << pp.get_volume() << '\t'
					<< pp.get_mean_fractal_dimension() << '\t'
					<< dt << '\t'
                    << clock.interval() / PRINT_EVERY << std::endl;

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
				
				clock.start();
            }

			// Save particles sizes for PSD
			//if (save_step>=PSD_DATA) {

			//	clock.stop();

			//	std::valarray<double> particles_sizes = pp.get_particles_sizes();
			//	std::valarray<double> aggregates_sizes = pp.get_aggregates_sizes();

			//	// print particles sizes
			//	if (particles_sizes.size() > 0) {

			//		for (int i = 0; i < particles_sizes.size(); i++) {
			//			part_sizes_file << particles_sizes[i] << " ";
			//		}
			//		part_sizes_file << std::endl;
			//	}

			//	// print aggregates sizes
			//	if (aggregates_sizes.size() > 0) {
			//		for (int i = 0; i < aggregates_sizes.size(); i++) {
			//			agg_sizes_file << aggregates_sizes[i] << " ";
			//		}
			//		agg_sizes_file << std::endl;
			//	}

			//	save_step = 0.0;

			//	clock.start();
			//}

        } // end simulation main loop

		// close plot file
		plot_data.close();

		// Anyway at the end save aggregates and primary particles data for PSD
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

    }// end loop on streamlines

	global_psd_aggregates.close();
	global_psd_primary.close();

	system("PAUSE");

}

void particle_link_temp() {


	WallClock clock;

	// for testing purposes
	// Declare the path for the XML file
	//std::string xml_path = "streams/NanoDomeStreamlines.xml";
	std::string xml_path = "streams/NanoDomeStreamlines_Marcus_test.xml";
	//std::string xml_path = "streams/TestStream.xml";
	// Declare the path for the JSON file
	//std::string json_path = "streams/NanoDomeStreamlines.json";

	// Global PSD from all streamlines (primary particles)
	std::ofstream global_psd_primary;
	std::string global_psd_primary_file = "global_primary_PSD.dat";
	global_psd_primary.open(global_psd_primary_file);

	// Global PSD from all streamlines (aggregates)
	std::ofstream global_psd_aggregates;
	std::string global_psd_aggregates_file = "global_aggregates_PSD.dat";
	global_psd_aggregates.open(global_psd_aggregates_file);

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	// FOR NOW NO SPECIES MANAGEMENT FROM THE STREAMLINE
	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si, ar };


	/// Execute a simulation for each streamline in the XML file

	// Streamlines loop
	for (int s = 0; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;

		// create splines with samples
		l_streams_xml[s].create_splines();

		// streamline start time
		double start_time = l_streams_xml[s].get_minTime();
		// starting molar fraction value
		double start_mf = l_streams_xml[s].get_splineMF()->get_v(start_time);

		GasPhaseLink gp(l_streams_xml[s].get_splineT(),
			l_streams_xml[s].get_splineP(),
			gas_species, start_time, start_mf);

		// setup the nucleation theory we want to use
		ClassicalNucleationTheory cnt(si);

		// setup the particle phase we want to use
		PBMFractalParticlePhase<PBMAggregate<Particle>> pp(1.6, 1.0e-9);

		// set max and min number of aggregates
		pp.set_max_aggregates(500);
		pp.set_min_aggregates(490);

		// get the streamline starting time
		double t = start_time + 1.0e-10;

		// get the end time for the streamline
		double end_time = l_streams_xml[s].get_maxTime();

		int iter = 0;

		const int PRINT_EVERY = 10000;
		
		const int SAVE_EVERY = 2000;
		std::ofstream plot_data;
		std::string filename = "PBM_plot" + std::to_string(s) + "Streamline.dat";
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

		// PSD Data files
		const double PSD_DATA = 1.0e-5;
		double save_step = 0.0;
		// Prinary particles
		std::ofstream part_sizes_file;
		std::string part_size_filename = "particles_sizes" + std::to_string(s) + "Streamline.dat";
		part_sizes_file.open(part_size_filename);
		// Aggregates diameters
		std::ofstream agg_sizes_file;
		std::string agg_size_filename = "aggregates_sizes" + std::to_string(s) + "Streamline.dat";
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

			gp.timestep_temp_grad(t, dt / 2.0, { 0.0, 0.0 });
			pp.volume_expansion(dt / 2.0, gp);

			//t += dt;
			save_step += dt;
			iter++;

			if (counter_trigger(iter, PRINT_EVERY)) {

				clock.stop();

				double T = gp.get_T();
				double S = gp.get_S("Si");

				std::cout
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
					<< dt << '\t'
					<< clock.interval() / PRINT_EVERY << std::endl;

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

				clock.start();
			}

			// Save particles sizes for PSD
			if (save_step >= PSD_DATA) {

				clock.stop();

				std::valarray<double> particles_sizes = pp.get_particles_sizes();
				std::valarray<double> aggregates_sizes = pp.get_aggregates_sizes();

				// print particles sizes
				if (particles_sizes.size() > 0) {

					for (int i = 0; i < particles_sizes.size(); i++) {
						part_sizes_file << particles_sizes[i] << " ";
					}
					part_sizes_file << std::endl;
				}

				// print aggregates sizes
				if (aggregates_sizes.size() > 0) {
					for (int i = 0; i < aggregates_sizes.size(); i++) {
						agg_sizes_file << aggregates_sizes[i] << " ";
					}
					agg_sizes_file << std::endl;
				}

				save_step = 0.0;

				clock.start();
			}

		} // end simulation main loop

		  // close plot file
		plot_data.close();

		// Anyway at the end save aggregates and primary particles data for PSD
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

	}// end loop on streamlines

	// close global PSD files
	global_psd_aggregates.close();
	global_psd_primary.close();

	system("PAUSE");


}

/*---------- TEST PARALLEL IMPLEMENTATION -----------------------------*/
void particle_link_thd(streamline& s, int s_idx) {

    WallClock clock;

    // FOR NOW NO SPECIES MANAGEMENT FROM THE STREAMLINE
    Species si("Si");
    Species ar("Ar");

    // create the vector of species for the gas phase
    std::vector<Species> gas_species = { si, ar };

    // Create Spline for Temperature
    SplineLinear *sT = new SplineLinear(s.get_Time(), s.get_Temp());

    // Create Spline for Pressure
    SplineLinear *sP = new SplineLinear(s.get_Time(), s.get_Press());

    // Create list of Spline for species concentration
    //SplineLinear *sC_si = new SplineLinear(s.get_Time(), s.get_Molar());
    //std::vector<Spline*> l_sc; l_sc.push_back(sC_si);

    // Create list of Spline for species molar fraction
    SplineLinear *sC_f_si = new SplineLinear(s.get_Time(), s.get_Molar_Fraction());
    std::vector<Spline*> l_sc; l_sc.push_back(sC_f_si);

    // Dummy Spline for Ar
    SplineLinear *sC_ar = new SplineLinear();
    l_sc.push_back(sC_ar);

	double T_start = sT->get_v(sT->get_x_min());

    // set pressure, temperature, species and initial molar concentration
    GasPhaseLink gp(sT, sP, gas_species, l_sc, sT->get_x_min());

    // setup the nucleation theory we want to use
    ClassicalNucleationTheory cnt(si);

    // setup the particle phase we want to use
    PBMFractalParticlePhase<PBMAggregate<Particle>> pp(1.6, 1.0e-9);

	// set max and min number of aggregates
	pp.set_max_aggregates(2000);
	pp.set_min_aggregates(1990);

    //double t = 0.0;
    // get the streamline starting time
    double t = sT->get_x_min();
    int iter = 0;

	const int PRINT_EVERY = 1000;

	const int SAVE_EVERY = 500;
	std::ofstream plot_data;
	std::string filename = "RESULTS/PBM_plot" + std::to_string(s.get_ID()) + "Streamline.dat";
	plot_data.open(filename);

	// PSD Data files
	const double PSD_DATA = 1.0e-5;
	double save_step = 0.0;
	// Prinary particles
	std::ofstream part_sizes_file;
	std::string part_size_filename = "RESULTS/particles_sizes" + std::to_string(s.get_ID()) + "Streamline.dat";
	part_sizes_file.open(part_size_filename);
	// Aggregates diameters
	std::ofstream agg_sizes_file;
	std::string agg_size_filename = "RESULTS/aggregates_sizes" + std::to_string(s.get_ID()) + "Streamline.dat";
	agg_sizes_file.open(agg_size_filename);

	clock.start();

    std::cout << "STREAMLINE: " << s.get_ID() << " STARTED" << std::endl;

    clock.start();

    // loop over timesteps
    double end_time = sT->get_x_max();
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
		if (std::isnan(dt)) {
			std::cout << "dt!->rho:  " << rho << " total rate: " << R_tot << "-log(rho)" << -log(rho) << std::endl;
		}

		// Strang first step
		//gp.timestep(dt / 2.0, dTdt, 0);
		t += (dt / 2.0); // First half-step
		if (t > end_time)
			break;

		gp.timestep(t, dt / 2.0);
		pp.volume_expansion(dt / 2.0, gp);

		// Strang second step
		g_si += pp.timestep(dt, gp, cnt);
		//gp.timestep(dt, 0, 0, { -g_si,0.0 });
		gp.timestep(t, dt / 2.0);

		// Strang third step
		//gp.timestep(dt / 2.0, dTdt, 0);
		t += (dt / 2.0); // Second half-step
		if (t > end_time)
			break;

		gp.timestep(t, dt / 2.0);
		pp.volume_expansion(dt / 2.0, gp);

		//t += dt;
		save_step += dt;
		iter++;

		//if (counter_trigger(iter, PRINT_EVERY)) {

		//	clock.stop();

		//	double T = gp.get_T();
		//	double S = gp.get_S("Ti");

		//	std::cout
		//		<< t << '\t'
		//		<< gp.get_T() << '\t'
		//		<< gp.get_S("Ti") << '\t'
		//		<< cnt.nucleation_rate(T, S) << '\t'
		//		<< gp.get_n("Ti") << '\t'
		//		<< cnt.stable_cluster_size(T, S) << '\t'
		//		<< "AVG Part Num: " << pp.get_mean_particles_number() << '\t'
		//		<< "S_level: " << pp.get_mean_sintering_level() << '\t'
		//		<< pp.get_aggregates_mean_spherical_diameter() << '\t'
		//		<< pp.get_aggregates_number() << '\t'
		//		<< pp.get_aggregates_density() << '\t'
		//		<< pp.get_volume() << '\t'
		//		<< pp.get_mean_fractal_dimension() << '\t'
		//		<< dt << '\t'
		//		<< clock.interval() / PRINT_EVERY << std::endl;

		//	// Print file

		//	clock.start();
		//}

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

			clock.start();
		}

		// Save particles sizes for PSD
		if (save_step >= PSD_DATA) {

			clock.stop();

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

			save_step = 0.0;

			clock.start();
		}

	} // end simulation main loop

	  // close plot file
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

	std::cout << "STREAMLINE: " << s_idx << " ENDED" << std::endl;

}

void particle_parallel() {

    // READ STREAMLINES FROM THE XML FILE

    // for testing purposes
    // Declare the path for the XML file
    std::string xml_path = "streams/NanoDomeStreamlines.xml";
    // Declare the path for the JSON file
    //std::string json_path = "streams/NanoDomeStreamlines.json";

    /*---------------------- XML TESTING -----------------------------*/

    // Create the vector containing all the streamlines in the file
    std::vector<streamline> l_streams_xml;

    // Inizialize the data structure to manage the XML file
    streamlinefileXML s_xml(xml_path, l_streams_xml);

    // Create the streamline
    s_xml.read_Streamlines();
	
	// Parallel execution of the simulations
    WallClock clock;
    omp_set_num_threads(4);
    int s;

    clock.start();

#pragma omp parallel for private(s)
    for (s = 0; s < l_streams_xml.size(); s++) {

       particle_link_thd(l_streams_xml[s], s);
	}

    clock.stop();
    std::cout << "ELAPSED TIME: " << clock.interval() << std::endl;


    system("PAUSE");


}

/*---------- END TEST PARALLEL IMPLEMENTATION -----------------------------*/

// Langevin particle approach
#include "dynamicparticle.h"
#include "aggregate/rattleaggregate.h"
#include "particlephase/constrainedlangevinparticlephase.h"
#include <experimental/filesystem>

void langevinPP() {

    WallClock clock;
    int PRINT_STEP = 40;
	double SAVE_SNAPSHOT = 1.0e-8;
	std::string vtk_path = "E/vtk/";

    const int PRINT_EVERY = 1000;
    const int SAVE_EVERY = 100;

    Species si("Si");
    Species ar("Ar");

    // create the vector of species for the gas phase
    std::vector<Species> gas_species = { si, ar };

    // gas phase info
    double p = 1.01e5; // pressure [Pa]
    double dTdt = -1e7; // temperature gradient [K/s]

    // Simulation starting paramenters
    double T_start = 3000.0;
    double T_end = 700.0;
    double V_start = 9.0e-18;
    double dt = 0.1e-9;

    // set pressure, temperature, species and initial molar concentration
    GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75 });

    // setup the nucleation theory we want to use
    ClassicalNucleationTheory cnt(si);

    // setup the particle phase
    ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> pp(V_start);

    double elapsed_time = 0.0; // Simulation time
    int iterations = 0; // Iterations
	double snap_count = 0.0; // Counter for saving VTK file

    double T = T_start;

    // Simulation Main Cycle
    while (T > T_end) {

        // check the smallest particle in the system
        double d_min = pp.get_particles_smallest_diameter();

        // calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
        double dt_max_lang = d_min*si.get_bulk_density(gp.get_T())/gp.get_gas_flux();
        // calculate dt max to have v*dt < d/2 for the smallest particle
        double dt_max_coll = sqrt(M_PI*si.get_bulk_density(gp.get_T()) *pow(d_min,5)/(24*3*K_BOL*gp.get_T()));

        if(pp.get_aggregates_number()>=1) {
            dt = std::min(dt_max_coll,dt_max_lang);
            if(dt>1e-10) dt = 1.0e-10;
        }

        // particle phase timestep
        double g_si = pp.timestep(dt, gp, cnt, T);

        // gas phase time step
        gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

        // Update elapsed time and iterations
        elapsed_time += dt;
		snap_count += dt;
        iterations++;

        // Print snapshot
        if (counter_trigger(iterations,PRINT_EVERY)) {
            std::cout << elapsed_time << '\t'
                << dt_max_coll << '\t'
                << dt_max_lang << '\t'
                << "T[K]: " << gp.get_T() << '\t'
                << "V[m^3]"<< pp.get_volume() << '\t'
                //<< J_si << '\t'
                //<< ns << '\t'
                //<< j_si << '\t'
                << "|N|: " << pp.get_aggregates_number()
                //<< "|N|: " << pp.get_aggregates_cardinality()
                << std::endl;
		}

		// Save VTK
		if (snap_count >= SAVE_SNAPSHOT && pp.get_aggregates_number() > 0) {
			pp.save_vtk(iterations, vtk_path);
			snap_count = 0.0;
		}

    }

    system("PAUSE");

}

void langevinPP_link() {

	WallClock clock;
	int PRINT_STEP = 40;
	double SAVE_SNAPSHOT = 1.0e-8;

	// for testing purposes
	// Declare the path for the XML file
	std::string xml_path = "streams/NanoDomeStreamlines.xml";
    std::string root = "vtk/";

	// Declare the path for the JSON file
	//std::string json_path = "streams/NanoDomeStreamlines.json";

	/*---------------------- XML TESTING -----------------------------*/

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	// FOR NOW NO SPECIES MANAGEMENT FROM THE STREAMLINE
	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si, ar };


	/// Execute a simulation for each streamline in the read file

	for (int s = 0; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;
		// Creates Streamline folder
		std::string stream_folder = root + "Streamline" + std::to_string(s) + "/";
		
		// Create the folder for the streamline
		std::experimental::filesystem::create_directory(stream_folder);

		// Create splines interpolating the samples extracted
		l_streams_xml[s].create_splines();

		std::vector<Spline*> l_sc;
		l_sc.push_back(l_streams_xml[s].get_splineMF());

		// Dummy Spline for Ar
		SplineLinear *sC_ar = new SplineLinear();
		l_sc.push_back(sC_ar);

		// simulation timestep
		double dt = 1.0e-10;
		// Temperature
		double T = 0.0;

		GasPhaseLink gp(l_streams_xml[s].get_splineT(),
			l_streams_xml[s].get_splineP(),
			gas_species, l_sc,
			l_streams_xml[s].get_minTime());


		// setup the nucleation theory we want to use
		ClassicalNucleationTheory cnt(si);

		// Starting pressure 
		double V_start = 9.0e-18;

		// setup the particle phase
		ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> pp(V_start);

		double elapsed_time = l_streams_xml[s].get_minTime(); // Simulation time
		double end_time = l_streams_xml[s].get_maxTime(); // Simulation end Time
		int iterations = 0; // Iterations
		double snap_count = 0.0; // Counter for saving VTK file

		const int PRINT_EVERY = 100;

		const int SAVE_EVERY = 30;
		std::ofstream plot_data;
		std::string filename = "CG_plot" + std::to_string(s) + "Streamline.dat";
		plot_data.open(filename);

		// Simulation Main Cycle
		while (elapsed_time < end_time) {

			// check the smallest particle in the system
			double d_min = pp.get_particles_smallest_diameter();

			// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
			double dt_max_lang = d_min*si.get_bulk_density(gp.get_T()) / gp.get_gas_flux();
			// calculate dt max to have v*dt < d/2 for the smallest particle
			double dt_max_coll = sqrt(M_PI*si.get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

			if (pp.get_aggregates_number() >= 1) {
				dt = std::min(dt_max_coll, dt_max_lang);
				if (dt > 1e-10) dt = 1.0e-10;
			}

			T = gp.get_T();

			// particle phase timestep
			double g_si = pp.timestep(dt, gp, cnt, T);

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

				std::cout << elapsed_time << '\t'
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
					<< clock.interval() / PRINT_EVERY << std::endl;

				// Print file

				clock.start();
			}

			// Print file for plotting
			if (counter_trigger(iterations, SAVE_EVERY)) {

				clock.stop();

				double T = gp.get_T();
				double S = gp.get_S("Si");

				plot_data << elapsed_time << '\t'
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

	}

	system("PAUSE");
}

void langevinLC() {

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si, ar };

	// gas phase info
	double p = 1.01e5; // pressure [Pa]
	double dTdt = -1e7; // temperature gradient [K/s]

	// Simulation starting paramenters
	double T_start = 3000.0;
	double T_end = 300.0;
	double time_end = 0.002;
    double V_start = 1.0e-17;
	double dt = 0.1e-9;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75 });

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(si);

	// setup the particle phase
	ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> pp(V_start);
	// min and max number of particles
    pp.set_max_agg_number(100);
    pp.set_min_agg_number(90);

	double elapsed_time = 0.0; // Simulation time
	int iterations = 0; // Iterations
	double snap_count = 0.0; // Counter for saving VTK file

	double T = T_start;

	WallClock clock;
	double SAVE_SNAPSHOT = 5.0e-8;
	const int PRINT_EVERY = 1000;
	const int SAVE_EVERY = 500;

	std::ofstream plot_data;
	std::string filename = "Langevin_plot.dat";
    std::string vtk_path = "/mnt/raid1/strauba_data/vtk/";
	plot_data.open(filename);

	// PSD Data files
	const double PSD_DATA = 1.0e-5;
	double save_step = 0.0;
	// Prinary particles
	std::ofstream part_sizes_file;
	std::string part_size_filename = "particles_sizes.dat";
	part_sizes_file.open(part_size_filename);
	// Aggregates diameters
	std::ofstream agg_sizes_file;
	std::string agg_size_filename = "aggregates_sizes.dat";
	agg_sizes_file.open(agg_size_filename);

	// Simulation Main Cycle
	while (elapsed_time < time_end) {

		// check the smallest particle in the system
		double d_min = pp.get_particles_smallest_diameter();
		// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
		double dt_max_lang = d_min*si.get_bulk_density(gp.get_T()) / gp.get_gas_flux();
		// calculate dt max to have v*dt < d/2 for the smallest particle
		double dt_max_coll = sqrt(M_PI*si.get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

		if (pp.get_aggregates_number() >= 1) {
			dt = std::min(dt_max_coll, dt_max_lang);
			if (dt > 1.0e-9) 
				dt = 1.0e-9;
			if (dt < 1.0e-12) {
				std::cout << "Stability: " << dt_max_lang << " Collision: " << dt_max_coll << " Smallest Diameter: " << d_min << std::endl;
			}
		}

		// particle phase timestep
		double g_si = pp.timestep_lc(dt, gp, cnt, T);

		// gas phase time step
		if (T <= 300) dTdt = 0.0;
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
				<< elapsed_time << '\t'									// time
				<< gp.get_T() << '\t'							// temperature
				<< gp.get_S("Si") << '\t'						// supersaturation (S)
				<< cnt.nucleation_rate(T, S) << '\t'			// J
				<< gp.get_n("Si") << '\t'						// ns
				<< cnt.stable_cluster_size(T, S) << '\t'		// j
				<< pp.get_mean_particles_number() << '\t'		// N_m
				<< pp.get_mean_sintering_level() << '\t'		//
				<< pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< pp.get_aggregates_number() << '\t'
				<< pp.get_aggregates_density() << '\t'
				<< dt << '\t'
				<< pp.get_volume()<<'\t'
				<< pp.get_mean_fractal_dimension() << '\t'
				<< clock.interval() / PRINT_EVERY << std::endl;

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

		// Save particles sizes for PSD
		if (save_step >= PSD_DATA) {

			clock.stop();

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

			save_step = 0.0;

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

	system("PAUSE");
}

void langevinLC_link() {

	WallClock clock;
	
	// for testing purposes
	// Declare the path for the XML file
	std::string xml_path = "streams/NanoDomeStreamlines.xml";
	std::string root = "vtk/";
	// Declare the path for the JSON file
	//std::string json_path = "streams/NanoDomeStreamlines.json";

	/*---------------------- XML TESTING -----------------------------*/

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	// FOR NOW NO SPECIES MANAGEMENT FROM THE STREAMLINE
	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si, ar };


	/// Execute a simulation for each streamline in the read file

	for (int s = 0; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;
		
		// Creates Streamline folder
		std::string stream_folder = root + "Streamline" + std::to_string(s) + "/";

		std::experimental::filesystem::create_directory(stream_folder);

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
		ClassicalNucleationTheory cnt(si);

		// Starting pressure 
		double V_start = 5.0e-17;

		// setup the particle phase
		ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> pp(V_start);

		double elapsed_time = l_streams_xml[s].get_minTime(); // Simulation time
		double end_time = l_streams_xml[s].get_maxTime(); // Simulation end Time

		int iterations = 0; // Iterations
		double snap_count = 0.0; // Counter for saving VTK file

		const int PRINT_EVERY = 500;
		double SAVE_SNAPSHOT = 5.0e-8;
		const int SAVE_EVERY = 300;
		std::ofstream plot_data;
		std::string filename = "CG_plot" + std::to_string(s) + "Streamline.dat";
		plot_data.open(filename);

		// Simulation Main Cycle
        while (elapsed_time < end_time) {

			// check the smallest particle in the system
			double d_min = pp.get_particles_smallest_diameter();

			// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
			double dt_max_lang = d_min*si.get_bulk_density(gp.get_T()) / gp.get_gas_flux();
			// calculate dt max to have v*dt < d/2 for the smallest particle
			double dt_max_coll = sqrt(M_PI*si.get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

			if (pp.get_aggregates_number() >= 1) {
				dt = std::min(dt_max_coll, dt_max_lang);
				if (dt > 1e-8) dt = 1.0e-8;
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
					<< dt << '\t'
					<< pp.get_mean_fractal_dimension() << '\t'
					<< clock.interval() / PRINT_EVERY << std::endl;

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

	}

	system("PAUSE");


}

void PBM_Langevin() {

	WallClock clock;

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si, ar };

	// gas phase info
	double p = 1.01e5; // pressure [Pa]
	double dTdt = -1.0e7; // temperature gradient [K/s]

	// Temperature Data
	double T_start = 3000.0;
	double T_end = 300.0;
	double T;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75 });

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(si);

	// setup the particle phase we want to use
	PBMFractalParticlePhase<PBMAggregate<Particle>> pbm_pp(2.8, 1e-13);

	double t = 0.0;
	int iter = 0;

	const int PRINT_EVERY = 1000;
	const int SAVE_EVERY = 50;

	std::ofstream plot_data;
	std::string filename = "LangevinPBM_plot.dat";
	std::string vtk_path = "E:/vtk/";
	plot_data.open(filename);

	clock.start();

	// loop over timesteps
	while (pbm_pp.get_mean_particles_number() <= 1.1) {

		if (gp.get_T()<300) dTdt = 0;

		// species source term for the gas phase
		double g_si = 0.0;

		// calculate the timestep using an exponential waiting time
		double R_tot = pbm_pp.get_total_processes_rate(gp, cnt);
		double rho = ndm::uniform_double_distr(ndm::rand_gen);

		// exponential waiting time
		double dt = -log(rho) / R_tot;

		// Strang first step
		gp.timestep(dt / 2.0, dTdt, 0);
		pbm_pp.volume_expansion(dt / 2.0, gp);

		// Strang second step
		g_si += pbm_pp.timestep(dt, gp, cnt);
		gp.timestep(dt, 0, 0, { -g_si,0.0 });

		// Strang third step
		gp.timestep(dt / 2.0, dTdt, 0);
		pbm_pp.volume_expansion(dt / 2.0, gp);

		t += dt;
		iter++;

		if (counter_trigger(iter, PRINT_EVERY)) {

			clock.stop();

			double T = gp.get_T();
			double S = gp.get_S("Si");

			std::cout 
				<< t << '\t'
				<< gp.get_T() << '\t'
				<< gp.get_S("Si") << '\t'
				<< cnt.nucleation_rate(T, S) << '\t'
				<< gp.get_n("Si") << '\t'
				<< cnt.stable_cluster_size(T, S) << '\t'
				<< pbm_pp.get_mean_particles_number() << '\t'
				<< pbm_pp.get_mean_sintering_level() << '\t'
				<< pbm_pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< pbm_pp.get_aggregates_number() << '\t'
				<< pbm_pp.get_aggregates_density() << '\t'
				<< clock.interval() / PRINT_EVERY << std::endl;

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
				<< pbm_pp.get_mean_particles_number() << '\t'
				<< pbm_pp.get_mean_sintering_level() << '\t'
				<< pbm_pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< pbm_pp.get_aggregates_number() << '\t'
				<< pbm_pp.get_aggregates_density() << '\t'
				<< pbm_pp.get_volume() << '\t'
				<< pbm_pp.get_mean_fractal_dimension() << '\t'
				<< clock.interval() / PRINT_EVERY << std::endl;

			// Print file

			clock.start();
		}
	}

	// END PBM
	std::cout << "END PBM" << std::endl
		<< "Number of aggregates: " << pbm_pp.get_aggregates_number() << std::endl
		<< "Mean Sintering Level: " << pbm_pp.get_mean_sintering_level() << std::endl
		<< "Volume: " << pbm_pp.get_volume() << std::endl
		<< "Fractal dimension (3 = SPHERE): " << pbm_pp.get_mean_fractal_dimension() << std::endl
		<< "Mean Number of particles for aggregate: " << pbm_pp.get_mean_particles_number() << std::endl;

	// START LANGEVIN
	std::cout << "START LANGEVIN" << std::endl;


	T = gp.get_T();

	// setup the particle phase
	//ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> cg_pp(pbm_pp, pbm_pp.get_volume(), T);
	ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> cg_pp(pbm_pp, pbm_pp.get_volume(), T);

	double elapsed_time = t; // Simulation time
	double end_time = 0.01;
	int iterations = 0; // Iterations
	double snap_count = 0.0; // Counter for saving VTK file

	double SAVE_SNAPSHOT = 5.0e-8;

	

	// Simulation Main Cycle
	while (elapsed_time < end_time) {

		// check the smallest particle in the system
		double d_min = cg_pp.get_particles_smallest_diameter();

		// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
		double dt_max_lang = d_min*si.get_bulk_density(gp.get_T()) / gp.get_gas_flux();
		// calculate dt max to have v*dt < d/2 for the smallest particle
		double dt_max_coll = sqrt(M_PI*si.get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

		double dt = 0.0;
		if (cg_pp.get_aggregates_number() >= 1) {
			dt = std::min(dt_max_coll, dt_max_lang);
			if (dt > 5.0e-09) dt = 5.0e-09;
		}

		// particle phase timestep
		double g_si = cg_pp.timestep_lc(dt, gp, cnt, T);

		if (gp.get_T() < 300)
			dTdt = 0.0;

		// gas phase time step
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		// Update elapsed time and iterations
		elapsed_time += dt;
		snap_count += dt;
		iterations++;

		// Print snapshot
		if (iterations%PRINT_EVERY == 0 && cg_pp.get_aggregates_number() > 0) {
			clock.stop();

			double T = gp.get_T();
			double S = gp.get_S("Si");

			std::cout 
				<< elapsed_time << '\t'									// time
				<< gp.get_T() << '\t'							// temperature
				<< gp.get_S("Si") << '\t'						// supersaturation (S)
				<< cnt.nucleation_rate(T, S) << '\t'			// J
				<< gp.get_n("Si") << '\t'						// ns
				<< cnt.stable_cluster_size(T, S) << '\t'		// j
				<< cg_pp.get_mean_particles_number() << '\t'		// N_m
				<< cg_pp.get_mean_sintering_level() << '\t'		//
				<< cg_pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< cg_pp.get_aggregates_number() << '\t'
				<< cg_pp.get_aggregates_density() << '\t'
				<< dt << '\t'
				<< cg_pp.get_volume() << '\t'
				<< cg_pp.get_mean_fractal_dimension() << '\t'
				<< clock.interval() / PRINT_EVERY << std::endl;

			clock.start();
		}

		// Save VTK
		if (snap_count >= SAVE_SNAPSHOT && cg_pp.get_aggregates_number() > 0) {
			cg_pp.save_vtk(iterations, vtk_path);
			snap_count = 0.0;
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
				<< cg_pp.get_mean_particles_number() << '\t'
				<< cg_pp.get_mean_sintering_level() << '\t'
				<< cg_pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< cg_pp.get_aggregates_number() << '\t'
				<< cg_pp.get_aggregates_density() << '\t'
				<< cg_pp.get_volume() << '\t'
				<< cg_pp.get_mean_fractal_dimension() << '\t'
				<< clock.interval() / PRINT_EVERY << std::endl;


			clock.start();
		}

	}

	plot_data.close();

	/// Print fractal mean fractal dimension
	std::cout<<"MEAN FRACTAL DIMENSION: "<<cg_pp.get_mean_fractal_dimension()<<std::endl;

	system("PAUSE");


}

void PBM_Langevin_link() {

	WallClock clock;

	// for testing purposes
	// Declare the path for the XML file
	std::string xml_path = "streams/NanoDomeStreamlinesStream5.xml";
	std::string root = "E:/vtk/";
	// Declare the path for the JSON file
	//std::string json_path = "streams/NanoDomeStreamlines.json";

	/*---------------------- XML TESTING -----------------------------*/

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	// FOR NOW NO SPECIES MANAGEMENT FROM THE STREAMLINE
	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si, ar };


	/// Execute a simulation for each streamline in the read file

	for (int s = 0; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;

		// Creates Streamline folder
		std::string stream_folder = root + "Streamline" + std::to_string(s) + "/";

		// Create the directory for the streamline(Windows Only)

		//_mkdir(stream_folder.c_str());
		std::experimental::filesystem::create_directory(stream_folder);

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
		ClassicalNucleationTheory cnt(si);

		// setup the particle phase we want to use
		PBMFractalParticlePhase<PBMAggregate<Particle>> pbm_pp(1.6, 1.0e-9);

		// get the streamline starting time
		double t = l_streams_xml[s].get_minTime();

		// get the end time for the streamline
		double end_time = l_streams_xml[s].get_maxTime();

		int iter = 0;

		const int PRINT_EVERY = 1000;

		const int SAVE_EVERY = 500;
		std::ofstream plot_data;
		std::string filename = "PBM_Langevin_plot" + std::to_string(s) + "Streamline.dat";
		plot_data.open(filename);


		clock.start();
		
		// PBM MAIN CYCLE
		while (pbm_pp.get_mean_particles_number() <= 1.1) {

			//if (gp.get_T()<300) dTdt = 0;

			// species source term for the gas phase
			double g_si = 0.0;

			// calculate the timestep using an exponential waiting time
			double R_tot = pbm_pp.get_total_processes_rate(gp, cnt);
			double rho = ndm::uniform_double_distr(ndm::rand_gen);

			// exponential waiting time
			double dt = -log(rho) / R_tot;

			// Strang first step
			//gp.timestep(dt / 2.0, dTdt, 0);
			t += (dt / 2.0); // First half-step
			if (t > end_time)
				break;

			gp.timestep(t, dt / 2.0);
			pbm_pp.volume_expansion(dt / 2.0, gp);

			// Strang second step
			g_si += pbm_pp.timestep(dt, gp, cnt);
			//gp.timestep(dt, 0, 0, { -g_si,0.0 });
			gp.timestep(t, dt / 2.0);

			// Strang third step
			//gp.timestep(dt / 2.0, dTdt, 0);
			t += (dt / 2.0); // Second half-step
			if (t > end_time)
				break;

			gp.timestep(t, dt / 2.0);
			pbm_pp.volume_expansion(dt / 2.0, gp);

			//t += dt;
			iter++;

			if (counter_trigger(iter, PRINT_EVERY)) {

				clock.stop();

				double T = gp.get_T();
				double S = gp.get_S("Si");

				std::cout
					<< t << '\t'
					<< gp.get_T() << '\t'
					<< gp.get_S("Si") << '\t'
					<< cnt.nucleation_rate(T, S) << '\t'
					<< gp.get_n("Si") << '\t'
					<< cnt.stable_cluster_size(T, S) << '\t'
					<< pbm_pp.get_mean_particles_number() << '\t'
					<< pbm_pp.get_mean_sintering_level() << '\t'
					<< pbm_pp.get_aggregates_mean_spherical_diameter() << '\t'
					<< pbm_pp.get_aggregates_number() << '\t'
					<< pbm_pp.get_aggregates_density() << '\t'
					<< pbm_pp.get_volume() << '\t'
					<< pbm_pp.get_mean_fractal_dimension() << '\t'
					<< clock.interval() / PRINT_EVERY << std::endl;

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
					<< pbm_pp.get_mean_particles_number() << '\t'
					<< pbm_pp.get_mean_sintering_level() << '\t'
					<< pbm_pp.get_aggregates_mean_spherical_diameter() << '\t'
					<< pbm_pp.get_aggregates_number() << '\t'
					<< pbm_pp.get_aggregates_density() << '\t'
					<< pbm_pp.get_volume() << '\t'
					<< pbm_pp.get_mean_fractal_dimension() << '\t'
					<< clock.interval() / PRINT_EVERY << std::endl;

				// Print file

				clock.start();
			}

		}

		// END PBM
		std::cout << "END PBM" << std::endl
			<< "Number of aggregates: " << pbm_pp.get_aggregates_number() << std::endl
			<< "Mean Sintering Level: " << pbm_pp.get_mean_sintering_level() << std::endl
			<< "Volume: " << pbm_pp.get_volume() << std::endl
			<< "Fractal dimension " << pbm_pp.get_mean_fractal_dimension() << std::endl
			<< "Mean Number of particles for aggregate: " << pbm_pp.get_mean_particles_number() << std::endl;

		// START LANGEVIN
		std::cout << "START LANGEVIN" << std::endl;

		// Get the end Temperature of the PBM
		T = gp.get_T();

		// setup the Langevin particle phase
		ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> cg_pp(pbm_pp, pbm_pp.get_volume(), T);

		double elapsed_time = t; // Simulation time
		int iterations = 0; // Iterations
		double snap_count = 0.0; // Counter for saving VTK file

		double SAVE_SNAPSHOT = 5.0e-8;
		
		// Langevin MAIN CYCLE
		while (elapsed_time < end_time) {

			// check the smallest particle in the system
			double d_min = cg_pp.get_particles_smallest_diameter();

			// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
			double dt_max_lang = d_min*si.get_bulk_density(gp.get_T()) / gp.get_gas_flux();
			// calculate dt max to have v*dt < d/2 for the smallest particle
			double dt_max_coll = sqrt(M_PI*si.get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

			double dt = 0.0;
			if (cg_pp.get_aggregates_number() >= 1) {
				dt = std::min(dt_max_coll, dt_max_lang);
				if (dt > 5.0e-09) dt = 5.0e-09;
				if (dt < 1.0e-12) dt = 1.0e-12;
			}

			// particle phase timestep
			double g_si = cg_pp.timestep_lc(dt, gp, cnt, T);

			// gas phase time step
			gp.timestep(elapsed_time, dt);

			// Update elapsed time and iterations
			elapsed_time += dt;
			snap_count += dt;
			iterations++;

			// Print snapshot
			if (iterations%PRINT_EVERY == 0 && cg_pp.get_aggregates_number() > 0) {
				clock.stop();

				double T = gp.get_T();
				double S = gp.get_S("Si");

				std::cout
					<< elapsed_time << '\t'									// time
					<< gp.get_T() << '\t'							// temperature
					<< gp.get_S("Si") << '\t'						// supersaturation (S)
					<< cnt.nucleation_rate(T, S) << '\t'			// J
					<< gp.get_n("Si") << '\t'						// ns
					<< cnt.stable_cluster_size(T, S) << '\t'		// j
					<< cg_pp.get_mean_particles_number() << '\t'		// N_m
					<< cg_pp.get_mean_sintering_level() << '\t'		//
					<< cg_pp.get_aggregates_mean_spherical_diameter() << '\t'
					<< cg_pp.get_aggregates_number() << '\t'
					<< cg_pp.get_aggregates_density() << '\t'
					<< dt << '\t'
					<< cg_pp.get_volume() << '\t'
					<< cg_pp.get_mean_fractal_dimension() << '\t'
					<< clock.interval() / PRINT_EVERY << std::endl;

				clock.start();
			}

			// Save VTK
			if (snap_count >= SAVE_SNAPSHOT && cg_pp.get_aggregates_number() > 0) {
				cg_pp.save_vtk(iterations, stream_folder);
				snap_count = 0.0;
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
					<< cg_pp.get_mean_particles_number() << '\t'
					<< cg_pp.get_mean_sintering_level() << '\t'
					<< cg_pp.get_aggregates_mean_spherical_diameter() << '\t'
					<< cg_pp.get_aggregates_number() << '\t'
					<< cg_pp.get_aggregates_density() << '\t'
					<< cg_pp.get_volume() << '\t'
					<< cg_pp.get_mean_fractal_dimension() << '\t'
					<< clock.interval() / PRINT_EVERY << std::endl;


				clock.start();
			}

		}

		plot_data.close();

		/// Print fractal mean fractal dimension
		std::cout << "MEAN FRACTAL DIMENSION: " << cg_pp.get_mean_fractal_dimension() << std::endl;

	}
	system("PAUSE");
}

#include "cg_simulation_dataXML.h"
void test_hotstart() {

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si, ar };

	// gas phase info
	double p = 1.01e5; // pressure [Pa]
	double dTdt = -1e7; // temperature gradient [K/s]

						// Simulation starting paramenters
	double T_start = 3000.0;
	double T_end = 300.0;
	double T = 1200;
	double V_start = 5.0e-18;
	double dt = 0.1e-9;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T, gas_species, { 0.25, 0.75 });

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(si);

	// TEST LOAD FROM FILE
	CGSimulationDataXML simulation_Snap("test_read.vtp");

	auto agg_list = simulation_Snap.create_aggregates();

	ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> pp(agg_list, 5.0e-18, T);

	double elapsed_time = 0.0; // Simulation time
	int iterations = 0; // Iterations
	double snap_count = 0.0; // Counter for saving VTK file

	WallClock clock;
	double SAVE_SNAPSHOT = 5.0e-8;
	const int PRINT_EVERY = 500;
	const int SAVE_EVERY = 300;

	std::ofstream plot_data;
	std::string filename = "Langevin_plot.dat";
	std::string vtk_path = "E/vtk/";
	plot_data.open(filename);

	// Simulation Main Cycle
	while (T > T_end) {

		// check the smallest particle in the system
		double d_min = pp.get_particles_smallest_diameter();

		// calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
		double dt_max_lang = d_min*si.get_bulk_density(gp.get_T()) / gp.get_gas_flux();
		// calculate dt max to have v*dt < d/2 for the smallest particle
		double dt_max_coll = sqrt(M_PI*si.get_bulk_density(gp.get_T()) *pow(d_min, 5) / (24 * 3 * K_BOL*gp.get_T()));

		if (pp.get_aggregates_number() >= 1) {
			dt = std::min(dt_max_coll, dt_max_lang);
			if (dt > 1.0e-8) dt = 1.0e-8;
		}

		// particle phase timestep
		double g_si = pp.timestep_lc(dt, gp, cnt, T);

		// gas phase time step
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
				<< elapsed_time << '\t'									// time
				<< gp.get_T() << '\t'							// temperature
				<< gp.get_S("Si") << '\t'						// supersaturation (S)
				<< cnt.nucleation_rate(T, S) << '\t'			// J
				<< gp.get_n("Si") << '\t'						// ns
				<< cnt.stable_cluster_size(T, S) << '\t'		// j
				<< pp.get_mean_particles_number() << '\t'		// N_m
				<< pp.get_mean_sintering_level() << '\t'		//
				<< pp.get_aggregates_mean_spherical_diameter() << '\t'
				<< pp.get_aggregates_number() << '\t'
				<< pp.get_aggregates_density() << '\t'
				<< dt << '\t'
				<< pp.get_volume() << '\t'
				<< pp.get_mean_fractal_dimension() << '\t'
				<< clock.interval() / PRINT_EVERY << std::endl;

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

		// Save VTK
		if (snap_count >= SAVE_SNAPSHOT && pp.get_aggregates_number() > 0) {
			pp.save_vtk(iterations, vtk_path);
			snap_count = 0.0;
		}
	}
	// Close plot file
	plot_data.close();

	system("PAUSE");
}

void test_aggregate_copy() {

	Species s("Si");
	double n_mon = 100000000;

	std::valarray<double> x1 = { 11.0, 12.0, 13.0 };
	std::valarray<double> x2 = { 21.0, 22.0, 23.0 };
	std::valarray<double> x3 = { 31.0, 32.0, 33.0 };

	auto p1 = std::make_shared<DynamicParticle>(n_mon, s, x1);
	auto p2 = std::make_shared<DynamicParticle>(n_mon, s, x2);
	auto p3 = std::make_shared<DynamicParticle>(n_mon, s, x3);

	auto b1 = std::make_shared<ParticleBond<DynamicParticle>>(p1, p2);
	auto b2 = std::make_shared<ParticleBond<DynamicParticle>>(p2, p3);

	std::list<std::shared_ptr<DynamicParticle>> p_list;
	p_list.push_back(p1); p_list.push_back(p2); p_list.push_back(p3);

	std::list<std::shared_ptr<ParticleBond<DynamicParticle>>> b_list;
	b_list.push_back(b1); b_list.push_back(b2);

	auto p_agg = std::make_shared<RATTLEAggregate<DynamicParticle>>(p_list, b_list);

	std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> p_agg_list;

	p_agg_list.push_back(p_agg);

	auto p_agg2 = std::make_shared<RATTLEAggregate<DynamicParticle>>(*p_agg);
	p_agg_list.push_back(p_agg2);


}

#include "XML_config_file.h"
void test_xml_start(std::string xml_config_path) {

	XMLconfigfile Configuration(xml_config_path);
	Configuration.run_simulation();

	system("PAUSE");
}

#include "xml_templates.h"

void user_interface(int argc, char* argv[]) {

	std::string opt_1 = argv[1];
	if (opt_1 == "-config_file") {

		std::string config_file = argv[2];
		std::cout << "Configuration file: " << config_file << std::endl;
		// Start Simulation
		XMLconfigfile Configuration(config_file);
		Configuration.run_simulation();

		system("PAUSE");
	}
	else if (opt_1 == "-help") {
		std::cout << "NanoDome Simulation Environment, v.0.001 alpha" << std::endl;
		std::cout << "Available models:" << std::endl;
		std::cout << "Linked Simultions:" << std::endl;
		std::cout << '\t' << "Moments" << std::endl
			<< '\t' << "Population Based Model (PBM)" << std::endl
			<< '\t' << "Langevin" << std::endl
			<< '\t' << "Nodal Method (Testing)" << std::endl;
		std::cout << "Stand Alone Simultions:" << std::endl;
		std::cout << '\t' << "Moments" << std::endl
			<< '\t' << "Population Based Model (PBM)" << std::endl
			<< '\t' << "Langevin" << std::endl
			<< '\t' << "Nodal Method (Testing)" << std::endl;
		std::cout << "Template Creation: " << std::endl;
		std::cout << '\t' << "-template " << "[linked_moment || linked_pbm || linked_langevin || stand_alone_moment || stand_alone_pbm || stand_alone_langevin] " << "filename.xml" << std::endl;
		std::cout << "Execution from command line: executable_file -config_file filename.xml" << std::endl;


	}
	else if (opt_1 == "-template") {
		std::string opt_2 = argv[2];
		std::string filename = argv[3];

		if (opt_2 == "linked_moment") {
			linked_moments_template(filename);
		}
		else if (opt_2 == "linked_pbm") {
			linked_moments_template(filename);
		}
		else if (opt_2 == "linked_langevn") {
			linked_langevin_template(filename);
		}
		else if (opt_2 == "stand_alone_moment"){
			stand_alone_moment_template(filename);
		}
		else if (opt_2 == "stand_alone_pbm") {
			stand_alone_pbm_template(filename);
		}
		else if (opt_2 == "stand_alone_langevin") {
			stand_alone_langevin_template(filename);
		}
	}

	else {
		std::cout << "Option not recognized!" << std::endl;
	}



}

#include "sectional/sectional.h"
#include "dynamic_timestep.h"

void sectional() {

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si,ar };
	//std::vector<Species> gas_species = { al,ar };

	// simulation timestep
	const double dt = 1.0e-8;
	//const double dt = 1.0e-5;

	// gas phase info
	double p = 1.01e5; // pressure [Pa]
	//double dTdt = -1.0e6; // temperature gradient [K/s]
	double dTdt = 0.0; // temperature gradient [K/s]
	
	double T_start = 2800.0;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, { 0.0, 1.0});
	//GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75});
	//GasPhaseCV gp(p, T_start, gas_species, { 0.1, 0.9 });

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(si);
	
	// setup sectional method
	int bins = 60;
	double scale_f = 1.3;
	//double scale_f = 2.0;
	double starting_bin = 3.0*si.m_volume();
	//double starting_bin = 1.0e-27;
	SectionalModel sec_mod(bins, scale_f, starting_bin, si);

	// Create a dummy start distribution
	sec_mod.dummy_distribution(6.56e-26, 0.026, 5.17e19);

	double t = 0.;
	int iter = 0;
	WallClock clock;

	const int PRINT_EVERY = 1000;
	const int SAVE_EVERY = 1;
	std::ofstream plot_data;
	std::string filename = "SECTIONAL_plot.dat";

	// Clean BinValues
	std::string binval_path = "BinValues.dat";
	std::ofstream bin_val_f;
	bin_val_f.open(binval_path);
	bin_val_f.clear();
	bin_val_f.close();
	plot_data.open(filename);

	// Print headlines
	/*plot_data << "Time[sec]" << '\t'
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
*/
	sec_mod.print_bins_values();
	
	// loop over timesteps until the final temperature
	while (t < 0.001) {

		/*if (gp.get_T() >= 3500) dTdt = -dTdt;
		if (gp.get_T() <= T_start&& t > 0.0) break;*/
		//if (gp.get_T()>3300) dTdt = 0;
		//dTdt = Temp_gradient(t, 1);
		
		double T = gp.get_T();
		double ns = gp.get_n("Si");
		double S = ns / si.n_sat(T);
		
		//double S = ns / al.n_sat(T);
		double J = cnt.nucleation_rate(T, S);
		double j = cnt.stable_cluster_size(T, S);

		//if (T <= 2000) break;

		//dTdt = Temp_gradient(t, T, 1000);
		//dTdt = Temp_gradient(t, 2);

		// moment method timestep and species consumption retrieval
		double g_si = sec_mod.timestep(dt, gp, cnt, "Si");

		// updating the gas phase
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		++iter;
		t += dt;

		if (counter_trigger(iter, PRINT_EVERY)) {
			
			std::cout
				<< "time " << t << '\t'
				<< "Temp: " << T << '\t'
				<< "S: " << S << '\t'
				<< "J: " << J << '\t'
				<< "ns: " << ns << '\t'
				<< "j: " << j << '\t'
				<< "d: " << sec_mod.get_density() << '\t'
				<< "mean_d: " << sec_mod.get_mean_diameter() << '\t'
				<< "gi_cond " << sec_mod.get_gi_cond() << '\t'
				<< "gi_nucl " << sec_mod.get_gi_nucl() << '\t'
				<< std::endl;
		}

		// File savings
		if (counter_trigger(iter, SAVE_EVERY)) {

			sec_mod.print_bins(t);

			plot_data
				<< t << '\t'
				<< T << '\t'
				<< S << '\t'
				<< J << '\t'
				<< ns << '\t'
				<< j << '\t'
				<< 0.0 << '\t' //dummy particles number
				<< 0.0 << '\t' // dummy sintering level
				<< sec_mod.get_mean_diameter() << '\t'
				<< 0.0 << '\t' // dummy aggregates number
				<< sec_mod.get_density() << '\t'
				<< 0.0 << '\t' // dummy Volume
				<< 0.0 << '\t' // dummy fractal dimension
				<< 0.0 << '\t' // dummy interval time
				<< sec_mod.get_gi_cond()
				<< std::endl; 
		}
	}

	plot_data.close();

	sec_mod.print_bins(t);

	system("PAUSE");

}

void sectional_link() {

	// for testing purposes
	// Declare the path for the XML file
	//std::string xml_path = "streams/NanoDomeStreamlines_brux3.xml";
	std::string xml_path = "streams/NanoDomeStreamlines_test_2018_01_12.xml";
	//std::string xml_path = "streams/NanoDomeStreamlinesMarcus2.xml"; 
	// Declare the path for the JSON file
	//std::string json_path = "streams/NanoDomeStreamlines.json";

	/*---------------------- XML TESTING -----------------------------*/

	// Create the vector containing all the streamlines in the file
	std::vector<streamline> l_streams_xml;

	// Inizialize the data structure to manage the XML file
	streamlinefileXML s_xml(xml_path, l_streams_xml);

	// Create the streamline
	s_xml.read_Streamlines();

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si,ar };

	// Main cycle: a simulation for each streamline extracted from the CFD

	for (int s = 0; s < l_streams_xml.size(); s++) {

		std::cout << "STREAMLINE: " << s << std::endl;

		// Create splines interpolating the samples extracted
		l_streams_xml[s].create_splines();

		std::vector<Spline*> l_sc;
		l_sc.push_back(l_streams_xml[s].get_splineMF());

		// Dummy Spline for Ar
		SplineLinear *sC_ar = new SplineLinear();
		l_sc.push_back(sC_ar);

		// simulation timestep
		//double dt = 1e-7;
		double dt = 0.0;
		GasPhaseLink gp(l_streams_xml[s].get_splineT(),
			l_streams_xml[s].get_splineP(),
			gas_species, l_sc,
			l_streams_xml[s].get_minTime());


		// setup the nucleation theory we want to use
		ClassicalNucleationTheory cnt(si);

		// setup sectional method
		int bins = 60;
		double scale_f = 1.6;
		//double scale_f = 2.0;
		double starting_bin = 10.0*si.m_volume();
		//double starting_bin = 1.0e-27;
		SectionalModel sec_mod(bins, scale_f, starting_bin, si);

		// get starting time and end time
		double t = l_streams_xml[s].get_minTime() + dt;
		double end_t = l_streams_xml[s].get_maxTime();

		int iter = 0;
		WallClock clock;

		const int PRINT_EVERY = 1000;

		const int SAVE_EVERY = 1;
		std::ofstream plot_data;
		std::string filename = "SECTIONAL_plot_" + std::to_string(s) + "_Streamline.dat";

		// Clean BinValues
		std::string binval_path = "BinValues_"+ std::to_string(s) + "_Streamline.dat";
		std::ofstream bin_val_f;
		bin_val_f.open(binval_path);
		bin_val_f.clear();
		bin_val_f.close();
		plot_data.open(filename);

		// Print headlines
		/*plot_data << "Time[sec]" << '\t'
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
		*/
		// print first line of the bin values file with the volumes of the bin
		sec_mod.print_bins_values(s);

		/// Create dynamic timestep
		DynamicTimestep dyn_ts(l_streams_xml[s].get_Time(), l_streams_xml[s].get_splineT()->get_gradients());

		// loop over timesteps until the final temperature
		while (t < end_t) {

			double T = gp.get_T();
			double ns = gp.get_n("Si");
			double S = ns / si.n_sat(T);

			double J = cnt.nucleation_rate(T, S);
			double j = cnt.stable_cluster_size(T, S);

			dt = dyn_ts.get_t(t, dt);

			// moment method timestep and species consumption retrieval
			double g_si = sec_mod.timestep(dt, gp, cnt, "Si");
			
			// updating the gas phase
			gp.timestep(t, dt);

			++iter;
			t += dt;

			if (counter_trigger(iter, PRINT_EVERY)) {

				std::cout
					<< "time " << t << '\t'
					<< "Temp: " << T << '\t'
					<< "S: " << S << '\t'
					<< "J: " << J << '\t'
					<< "ns: " << ns << '\t'
					<< "j: " << j << '\t'
					<< "d: " << sec_mod.get_density() << '\t'
					<< "mean_d: " << sec_mod.get_mean_diameter() << '\t'
					<< "gi_cond " << sec_mod.get_gi_cond() << '\t'
					<< "gi_nucl " << sec_mod.get_gi_nucl() << '\t'
					<< "dt " << dt<< '\t'
					<< std::endl;
			}

			// File savings
			if (counter_trigger(iter, SAVE_EVERY)) {

				sec_mod.print_bins(t, s);

				plot_data
					<< t << '\t'
					<< T << '\t'
					<< S << '\t'
					<< J << '\t'
					<< ns << '\t'
					<< j << '\t'
					<< 0.0 << '\t' //dummy particles number
					<< 0.0 << '\t' // dummy sintering level
					<< sec_mod.get_mean_diameter() << '\t'
					<< 0.0 << '\t' // dummy aggregates number
					<< sec_mod.get_density() << '\t'
					<< 0.0 << '\t' // dummy Volume
					<< 0.0 << '\t' // dummy fractal dimension
					<< 0.0 << '\t' // dummy interval time
					<< sec_mod.get_gi_cond()
					<< "dt " << dt << '\t'
					<< std::endl;
			}
		}

		plot_data.close();

		sec_mod.print_bins(t, s);
	}

	system("PAUSE");


}

void sectional_Al() {

	Species al("Al");
	Species ar("Ar");

	// create the vector of species for the gas phase
	//std::vector<Species> gas_species = { si,ar };
	std::vector<Species> gas_species = { al,ar };

	// simulation timestep
	const double dt = 5.0e-8;
	//const double dt = 1.0e-5;

	// gas phase info
	double p = 1.01e5; // pressure [Pa]

	double dTdt;
	// temperature gradient [K/s] (Cooling)
	//double dTdt = -1e7; 
	//dTdt = -1.0e3;

	// temperature gradient [K/s] (Heating)
	//double dTdt = -1e7; 
	dTdt = +1.0e5;
	
	double T_start = 1773.0;

	// set pressure, temperature, species and initial molar concentration
	//GasPhaseCV gp(p, T_start, gas_species, { 1.75e-3/3.0, 1.0 - 1.75e-3 / 3.0 });
	GasPhaseCV gp(p, T_start, gas_species, { 1.75e-3, 1.0 - 1.75e-3 });
	//GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75});

	// setup the nucleation theory we want to use
	//ClassicalNucleationTheory cnt(si);
	ClassicalNucleationTheory cnt(al);

	// setup sectional method
	int bins = 40;
	//double scale_f = 1.8;
	double scale_f = 2.0;
	double starting_bin = 10.0*al.m_volume();
	//double starting_bin = 1.0e-27;
	//double starting_bin = 1.4137e-29;
	//SectionalModel sec_mod(bins, scale_f, starting_bin, si);
	SectionalModel sec_mod(bins, scale_f, starting_bin, al);

	double t = 0.;
	int iter = 0;

	const int PRINT_EVERY = 100;

	const int SAVE_EVERY = 100;
	std::ofstream plot_data;
	std::string filename = "SECTIONAL_plot_AL.dat";
	plot_data.open(filename);

	// Clean BinValues
	std::string binval_path = "BinValues.dat";
	std::ofstream bin_val_f;
	bin_val_f.open(binval_path);
	bin_val_f.clear();
	bin_val_f.close();

	// Print headlines
	/*plot_data << "Time[sec]" << '\t'
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
	*/
	sec_mod.print_bins_values();

	// loop over timesteps until the final temperature
	while (t < 0.5) {

		if (gp.get_T()<300) dTdt = 0;
		//if (gp.get_T()<300) break;
		if (gp.get_T() > 4000) break;


		double T = gp.get_T();
		//double ns = gp.get_n("Si");
		double ns = gp.get_n("Al");
		double S = ns / al.n_sat(T);

		//double S = ns / al.n_sat(T);
		double J = cnt.nucleation_rate(T, S);
		double j = cnt.stable_cluster_size(T, S);

		// moment method timestep and species consumption retrieval
		double g_si = sec_mod.timestep(dt, gp, cnt, "Al");

		// updating the gas phase
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		++iter;
		t += dt;

		if (counter_trigger(iter, PRINT_EVERY)) {
			
			std::cout
				<< "time " << t << '\t'
				<< "Temp: " << T << '\t'
				<< "S: " << S << '\t'
				<< "J: " << J << '\t'
				<< "ns: " << ns << '\t'
				<< "j: " << j << '\t'
				<< "d: " << sec_mod.get_density() << '\t'
				<< "mean_d: " << sec_mod.get_mean_diameter() << '\t'
				<< "gi_cond " << sec_mod.get_gi_cond() << '\t'
				<< "gi_nucl " << sec_mod.get_gi_nucl() << '\t'
				<< std::endl;
		}

		// File savings
		if (counter_trigger(iter, SAVE_EVERY)) {

			sec_mod.print_bins(t);

			plot_data
				<< t << '\t'
				<< T << '\t'
				<< S << '\t'
				<< J << '\t'
				<< ns << '\t'
				<< j << '\t'
				<< 0.0 << '\t' //dummy particles number
				<< 0.0 << '\t' // dummy sintering level
				<< sec_mod.get_mean_diameter() << '\t'
				<< 0.0 << '\t' // dummy aggregates number
				<< sec_mod.get_density() << '\t'
				<< 0.0 << '\t' // dummy Volume
				<< 0.0 << '\t' // dummy fractal dimension
				<< 0.0 << '\t' // dummy interval time
				<< sec_mod.get_gi_cond()
				<< std::endl;
		}
	}

	plot_data.close();

	sec_mod.print_bins(t);

	system("PAUSE");

}

void moment_Al() {

	Species al("Al");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { al,ar };

	// simulation timestep
	const double dt = 5.0e-6;

	// gas phase info
	double p = 1.01e5; // pressure [Pa]
	//double dTdt = -1e7; // temperature gradient [K/s]
	double dTdt = -1000.0;					//

	double T_start = 1772.0;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, { 1.75e-3/3.0, 1.0 - 1.75e-3/3.0});

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(al);

	// set the moment method with the condensing species
	//MomentModelFriedlander mm(si);
	MomentModelPratsinis mm(al);

	double t = 0.;
	int iter = 0;

	const int PRINT_EVERY = 100;

	const int SAVE_EVERY = 100;
	std::ofstream plot_data;
	std::string filename = "MOMENTS_plot_Al.dat";
	plot_data.open(filename);
	// Print headlines
	/*plot_data << "Time[sec]" << '\t'
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
	<< std::endl;*/



	// loop over timesteps until the final temperature
	while (t < 1.5) {

		if (gp.get_T()<300) dTdt = 0;

		double T = gp.get_T();
		double ns = gp.get_n("Al");
		double S = ns / al.n_sat(T);
		double J = cnt.nucleation_rate(T, S);
		double j = cnt.stable_cluster_size(T, S);

		if (S < 1.001 && S > 1.0)
			int a = 0;

		// moment method timestep and species consumption retrieval
		double g_si = mm.timestep(dt, gp, cnt);

		// updating the gas phase
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		++iter;
		t += dt;

		if (counter_trigger(iter, PRINT_EVERY)) {
			//gp.print();

			std::cout
				<< t << '\t'
				<< T << '\t'
				<< S << '\t'
				<< J << '\t'
				<< ns << '\t'
				<< j << '\t'
				<< mm.get_density() << '\t'
				<< mm.get_mean_diameter() << '\t'
				<< mm.get_cond_term()
				<< std::endl;
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
				<< 0.0 << '\t' // dummy interval time
				<< mm.get_cond_term()
				<< std::endl;
		}
	}

	plot_data.close();

	system("PAUSE");
}

void sectional_evap_SI() {

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si,ar };

	// simulation timestep
	const double dt = 1.0e-7;
	//const double dt = 1.0e-5;

	// gas phase info
	double p = 1.01e5; // pressure [Pa]
	//double dTdt = -1e7; // temperature gradient [K/s]
	double dTdt = 0.0;

	double T_start = 3000.0;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75 });
	//GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75});

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(si);

	// setup sectional method
	int bins = 40;
	double scale_f = 1.6;
	//double scale_f = 2.0;
	double starting_bin = 10.0* si.m_volume();
	//double starting_bin = 1.0e-27;
	SectionalModel sec_mod(bins, scale_f, starting_bin, si);


	double t = 0.;
	int iter = 0;

	const int PRINT_EVERY = 100;

	const int SAVE_EVERY = 100;
	std::ofstream plot_data;
	std::string filename = "SECTIONAL_plot_Evap_SI.dat";
	plot_data.open(filename);

	// Clean BinValues
	std::string binval_path = "BinValues.dat";
	std::ofstream bin_val_f;
	bin_val_f.open(binval_path);
	bin_val_f.clear();
	bin_val_f.close();

	// Print headlines
	/*plot_data << "Time[sec]" << '\t'
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
	*/
	sec_mod.print_bins_values();

	// loop over timesteps until the final temperature
	while (t < 0.005) {

		double T = gp.get_T();
		double ns = gp.get_n("Si");
		double S = ns / si.n_sat(T);

		//double S = ns / al.n_sat(T);
		double J = cnt.nucleation_rate(T, S);
		double j = cnt.stable_cluster_size(T, S);

		// moment method timestep and species consumption retrieval
		double g_si = sec_mod.timestep(dt, gp, cnt, "Si");

		// Sinusoidal Gradient
		dTdt = -2700.0*0.5 * 10000 * M_PI*0.5*std::sin(t * 10000 * M_PI*0.5);

		// updating the gas phase
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		++iter;
		t += dt;

		if (counter_trigger(iter, PRINT_EVERY)) {



			std::cout
				<< "time " << t << '\t'
				<< "Temp: " << T << '\t'
				<< "S: " << S << '\t'
				<< "J: " << J << '\t'
				<< "ns: " << ns << '\t'
				<< "j: " << j << '\t'
				<< "d: " << sec_mod.get_density() << '\t'
				<< "mean_d: " << sec_mod.get_mean_diameter() << '\t'
				<< "gi_cond " << sec_mod.get_gi_cond() << '\t'
				<< "gi_nucl " << sec_mod.get_gi_nucl() << '\t'
				<< std::endl;
		}

		// File savings
		if (counter_trigger(iter, SAVE_EVERY)) {

			sec_mod.print_bins(t);

			plot_data
				<< t << '\t'
				<< T << '\t'
				<< S << '\t'
				<< J << '\t'
				<< ns << '\t'
				<< j << '\t'
				<< 0.0 << '\t' //dummy particles number
				<< 0.0 << '\t' // dummy sintering level
				<< sec_mod.get_mean_diameter() << '\t'
				<< 0.0 << '\t' // dummy aggregates number
				<< sec_mod.get_density() << '\t'
				<< 0.0 << '\t' // dummy Volume
				<< 0.0 << '\t' // dummy fractal dimension
				<< 0.0 << '\t' // dummy interval time
				<< sec_mod.get_gi_cond()
				<< std::endl;
		}
	}

	plot_data.close();

	sec_mod.print_bins(t);

	system("PAUSE");

}

void moment_evap_SI() {

	Species si("Si");
	Species ar("Ar");

	// create the vector of species for the gas phase
	std::vector<Species> gas_species = { si,ar };

	// simulation timestep
	const double dt = 1.0e-8;

	// gas phase info
	double p = 1.01e5; // pressure [Pa]
	//double dTdt = -1e7; // temperature gradient [K/s]
	double dTdt = 0.0;

	// Starting Temperature
	double T_start = 3000.0;

	// set pressure, temperature, species and initial molar concentration
	GasPhaseCV gp(p, T_start, gas_species, { 0.25, 0.75 });

	// setup the nucleation theory we want to use
	ClassicalNucleationTheory cnt(si);

	// set the moment method with the condensing species
	//MomentModelFriedlander mm(si);
	MomentModelPratsinis mm(si);

	double t = 0.;
	int iter = 0;

	const int PRINT_EVERY = 1000;

	const int SAVE_EVERY = 100;
	std::ofstream plot_data;
	std::string filename = "MOMENTS_plot_evap_SI.dat";
	plot_data.open(filename);
	// Print headlines
	/*plot_data << "Time[sec]" << '\t'
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
	<< std::endl;*/



	// loop over timesteps until the final temperature
	while (t < 0.005) {

		double T = gp.get_T();
		double ns = gp.get_n("Si");
		double S = ns / si.n_sat(T);
		double J = cnt.nucleation_rate(T, S);
		double j = cnt.stable_cluster_size(T, S);

		// moment method timestep and species consumption retrieval
		double g_si = mm.timestep(dt, gp, cnt);

		// Sinusoidal temperature gradient
		dTdt = dTdt = -2700.0*0.5 * 10000 * M_PI*0.5*std::sin(t * 10000 * M_PI*0.5);

		// updating the gas phase
		gp.timestep(dt, dTdt, 0, { -g_si, 0.0 });

		++iter;
		t += dt;

		if (counter_trigger(iter, PRINT_EVERY)) {
			//gp.print();

			std::cout
				<< t << '\t'
				<< T << '\t'
				<< S << '\t'
				<< J << '\t'
				<< ns << '\t'
				<< j << '\t'
				<< mm.get_density() << '\t'
				<< mm.get_mean_diameter() << '\t'
				<< mm.get_cond_term()
				<< std::endl;
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
				<< 0.0 << '\t' // dummy interval time
				<< mm.get_cond_term()
				<< std::endl;
		}
	}

	plot_data.close();

	system("PAUSE");

}

int main(int argc, char* argv[]) {

	// Testing mode
	if (argc == 1) {
		moment();
		//particle();
		//langevinLC();
		//sectional();

		// Aluminum Tests
		//sectional_Al();
		//moment_Al();

		// Evaporation Tests
		//sectional_evap_SI();
		//moment_evap_SI();

		//moment_link();
		//particle_link();
		//langevinLC_link();
		//sectional_link();

		//moment_link_temp();
		//particle_link_temp();

		//langevinPP();
		//langevinPP_link();

		//load_test();
		//test_hotstart();
		//test_aggregate_copy();
		//test_xml_start("configuration_test_Link.xml");
		//test_xml_start("configuration_test_Link.xml");
		//test_xml_start("NanoDome_configuration.xml");


		//PBM_Langevin();
		//PBM_Langevin_link();

		/*moment_parallel();

		particle_parallel();*/
		//particle_parallel();

		//PBM_link_Zurich();
		//moments_link_zurich();
		//zurich_parallel();

	}
	else {
		// User interface mode
		user_interface(argc, argv);
		
	}

	return 0;
}
