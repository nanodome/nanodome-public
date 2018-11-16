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

#include "gasphaselink.h"


GasPhaseLink::GasPhaseLink(Spline* _T, Spline* _P, std::vector<Species> _species,
                           std::vector<Spline*> _C, double s_time): sT(_T), sP(_P), sC(_C){
    gamma = 1.0;

    T = sT->get_v(s_time);
    p = sP->get_v(s_time);
    species = _species;

    for (size_t i = 0; i<species.size(); ++i)
        hash[species[i].get_formula()] = i;

	c = { sC[0]->get_v(sC[0]->get_x_min()), 1- sC[0]->get_v(sC[0]->get_x_min()) };

    actual_time = s_time;

}

GasPhaseLink::GasPhaseLink(Spline* _T, Spline* _P, std::vector<Species> _species, double s_time, double nucl_species_start_m_frac)
	: sT(_T), sP(_P){

	gamma = 1.0;

	T = sT->get_v(s_time);
	p = sP->get_v(s_time);
	species = _species;

	for (size_t i = 0; i<species.size(); ++i)
		hash[species[i].get_formula()] = i;

	c = {nucl_species_start_m_frac, 1.0 - nucl_species_start_m_frac};

	actual_time = s_time;

}


void GasPhaseLink::timestep(double t, double dt) {

	double dTdt, dpdt; // T gradient and p gradient values

	// flag indicating if the next step is a significant step for the gradient
	bool sign_step = false;
	// getting the smallest value after t-dt
	double next_step = std::nextafter((t - dt), t);
	// if t (the next step in the simulation) is equals or bigger than next step, compute T and p gradients
	if (t - next_step > 0.0)
		sign_step = true;

	// Compute gradients between two values
	if (t > 0.0 && sign_step) { // the 
		// Temperature gradient
		double Tnum = sT->get_v(t) - sT->get_v(t - dt);
		double Tden = t - (t - dt);
		dTdt = Tnum / Tden;

		// Pressure gradient
		double pnum = sP->get_v(t) - sP->get_v(t - dt);
		double pden = t - (t - dt);
		dpdt = pnum / pden;
	}
	else // else set gradients to 0.0
	{
		std::cout << "Small timestep, gradients to 0.0" << std::endl;
		dTdt = dpdt = 0.0;
	}

    gamma = dTdt / T - dpdt / p;

    // gas phase pressure and temperature update
    T = sT->get_v(t);
    p = sP->get_v(t);
	MF = sC[0]->get_v(t);
	
	// Update molar fraction
	c = { sC[0]->get_v(t), 1 - sC[0]->get_v(t) };
}


void GasPhaseLink::timestep_temp_grad(double t, double dt, std::valarray<double> w) {

	double dTdt, dpdt; // T gradient and p gradient values

					   // flag indicating if the next step is a significant step for the gradient
	bool sign_step = false;
	// getting the smallest value after t-dt
	double next_step = std::nextafter((t - dt), t);
	// if t (the next step in the simulation) is equals or bigger than next step, compute T and p gradients
	if (t - next_step > 0.0)
		sign_step = true;

	// Compute gradients between two values
	if (t > 0.0 && sign_step) { // the 
								// Temperature gradient
		double Tnum = sT->get_v(t) - sT->get_v(t - dt);
		double Tden = t - (t - dt);
		dTdt = Tnum / Tden;

		// Pressure gradient
		double pnum = sP->get_v(t) - sP->get_v(t - dt);
		double pden = t - (t - dt);
		dpdt = pnum / pden;
	}
	else // else set gradients to 0.0
	{
		std::cout << "Small timestep, gradients to 0.0" << std::endl;
		dTdt = dpdt = 0.0;
	}

	double n = p / (K_BOL*T);
	double wtot = w.sum();

	std::valarray<double> ns = c*n;

	gamma = wtot / n + dTdt / T - dpdt / p;

	// gas phase pressure and temperature update
	T = sT->get_v(t);
	p = sP->get_v(t);

	// simple explicit ODE timestep
	// equation is solved for the number density
	for (std::size_t i = 0; i<c.size(); ++i)
		ns[i] += (w[i] - ns[i] * gamma) * dt;

	// molar concentration update
	n = p / (K_BOL*T);
	c = ns / n;


}

