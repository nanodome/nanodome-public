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

#ifndef DYNAMIC_TIMESTEP_H
#define DYNAMIC_TIMESTEP_H

#include <vector>

class DynamicTimestep {

	/// Timesteps for each gradient [sec.]
	std::vector<double> ts;

	/// Actual node
	int node_idx;

	/// Number of nodes
	int node_num;

	/// Gradients
	std::vector<double> gradients;

	/// time intervals [sec.]
	std::vector<double> t_nodes;

	/// Max Timestep
	double max_ts = 1.0e-6;

	/// Min Timestep
	double min_ts = 1.0e-9;

	/// Rule 1 for creating the timestep: ts = 1.0 / grad*10;
	double rule1(int _i);

public:

	/// Default Constructor
	///	\param: std::vector<double> _t_nodes time samples
	///	\param: std::vector<double> _grad gradients among time samples
	DynamicTimestep(std::vector<double> _t_nodes, std::vector<double> _grad);

	/// Returns the simulation timestep 
	///	\param double _t actual simulation time
	///	\param double _dt actual simulation timestep;
	double get_t(double _t, double _dt);

};



#endif

