#include "network_simulation.h"

NetworkSimulation::NetworkSimulation() {}

NetworkSimulation::NetworkSimulation(int _ID, double _T, double _p,
	std::list<std::string> _c_species, std::list<double> _c_species_m_fraction,
	std::list<std::string> _b_species, std::list<double> _b_species_m_fraction):
	ID(_ID), pressure(_p), Temperature(_T),
	c_species(_c_species), c_species_m_fraction(_c_species_m_fraction),
	b_species(_b_species), b_species_m_fraction(_b_species_m_fraction){}

void NetworkSimulation::send_particles() {

	auto b_it = send_buffer.begin();
	for (auto n_sim : neighbourgs) {
		n_sim->get_particles((*b_it));
		b_it++;
	}

}

void NetworkSimulation::get_particles(std::list<std::shared_ptr<Particle>> _part) {

	/// insert the received particles inside the receive buffer
	for (auto p : _part) {
		receive_buffer.push_back(p);
	}
}