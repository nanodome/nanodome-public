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

#include "dynamic_timestep.h"
#include <cmath>

DynamicTimestep::DynamicTimestep(std::vector<double> _t_nodes, std::vector<double> _grad): 
	t_nodes(_t_nodes), gradients(_grad){

	node_idx = 0;
	node_num = _t_nodes.size();

	ts.resize(node_num - 1);

	// Create timestep array
	for (int i = 0; i < node_num - 1; i++) {
		ts[i] = rule1(i);
	}
}

double DynamicTimestep::get_t(double _t, double _dt) {

	double dt = 0.0;



	if (_t + _dt > t_nodes[node_idx + 1] && node_idx < node_num) {
		dt = t_nodes[node_idx + 1] - _t;
		if (dt < min_ts)
			dt = min_ts;
		node_idx++;
		return dt;
	}
	else {
		return ts[node_idx];

	}

}

double DynamicTimestep::rule1(int _i) {

	double ts = 1.0 / (std::abs(gradients[_i]) * 10.0);

	if (ts < min_ts)
		return min_ts;
	if (ts > max_ts)
		return max_ts;

	return ts;

}
