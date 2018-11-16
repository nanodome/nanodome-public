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

#include "rattleaggregate.h"

template<typename P>
void RATTLEAggregate<P>::langevin_x_update(double dt, const GasPhase& gp, double T) {

    // for all DynamicParticles composing the aggregate, update x
    for (auto part: this->particles) {

        double dp = part->get_diameter();

        // free molecular friction coefficient [1/s] (in the reference article alpha = mass*gamma)
        double gamma = (4./3.) * M_PI * gp.get_gas_flux() * dp*dp / (4.*part->get_mass());

        part->langevin_verlet_x(dt, gamma, T);
    }
}


template<typename P>
void RATTLEAggregate<P>::langevin_v_update(double dt) {

    // for all DynamicParticles composing the aggregate, update v
    for (auto part : this->particles) {
        part->langevin_verlet_v(dt);
    }

}

template<typename P>
void RATTLEAggregate<P>::v_initialization(double T) {

    for (auto p = this->particles.begin(); p != this->particles.end(); p++) {
		(*p)->init_maxwellian_v(T);
	}

}

template<typename P>
void RATTLEAggregate<P>::langevin_F_update(double dt) {
    // TODO
}


template<typename P>
void RATTLEAggregate<P>::wall_reflection(double cube_side) {

    for (auto part : this->particles)
        part->wall_reflection(cube_side);
}

template<typename P>
void RATTLEAggregate<P>::periodic_domain(double _v_side_lenght) {

	for (auto part : this->particles)
		part->periodic_domain(_v_side_lenght);

}

template<typename P>
void RATTLEAggregate<P>::shift_coordinates(std::valarray<double> _shift) {

	for (auto p = this->particles.begin(); p != this->particles.end(); p++) {

		std::valarray<double> x = (*p)->get_x();
		x += _shift;
		(*p)->set_x(x);
	}

}

template<typename P>
void RATTLEAggregate<P>::change_center_of_mass(std::valarray<double> n_c_mass) {

	std::valarray<double> old_com = this->get_center_of_mass();
	std::valarray<double> diff = n_c_mass - old_com;

	for (auto p = this->particles.begin(); p != this->particles.end(); p++) {

		std::valarray<double> x = (*p)->get_x();
		x = x + diff;
		(*p)->set_x(x);
	}

}

template<typename P>
void RATTLEAggregate<P>::scale_coordinates(std::valarray<double> _scale) {

	for (auto p = this->particles.begin(); p != this->particles.end(); p++) {

		std::valarray<double> x = (*p)->get_x();
		x *= _scale;
		(*p)->set_x(x);
	}

}

template<typename P>
void RATTLEAggregate<P>::update_graph() {

#ifdef VERBOSE
	std::cout << "UPDATING GRAPH STRUCTURE" << std::endl;
#endif

    // Reset graph data
    graph.reset_graph();

    // Create new data structures
    graph.add_particles(this->particles, this->bonds);

    // Compute constraints
    graph.create_constraints();
}


template<typename P>
void RATTLEAggregate<P>::update_particles_cardinality() {
    actual_particles = this->particles.size();
}


template<typename P>
void RATTLEAggregate<P>::erase_graph() { graph.reset_graph(); }


template<typename P>
bool RATTLEAggregate<P>::SHAKE() {

	if (graph.SHAKE()) { 
		return true; 
	}
	else { 
		std::cout << "SHAKE for Aggregate: " << this->get_id()<<" Diverged, ELIMINATED!!"<<std::endl; 
		return false;
	}
}


template<typename P>
bool RATTLEAggregate<P>::check_coalescence() {
    if(actual_particles > this->particles.size() && this->particles.size()>0) {

#ifdef VERBOSE
        std::cout << "AGGREGATE ID: " << this->get_id() << std::endl;
        std::cout << "Actual Particles: " << actual_particles << std::endl;
        std::cout << "Listed Particles: " << this->particles.size() << std::endl;
#endif

        return true;

    } else {

        return false;
    }
}


template<typename P>
void RATTLEAggregate<P>::update_constraints() {
    graph.update_bonds(this->bonds);
}


template<typename P>
std::list<std::shared_ptr<P>> RATTLEAggregate<P>::get_particles()const { return this->particles; }


template<typename P>
void RATTLEAggregate<P>::print_test() {

    // print particle position
    for (auto it = this->particles.begin(); it != this->particles.end(); ++it) {
        (*it)->print_point();
    }
    // print vertices position of the graph structure
    //this->graph.print();
}

