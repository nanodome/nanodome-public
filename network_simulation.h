#ifndef NETWORK_SIMULATION_H
#define NETWORK_SIMULATION_H

#include "meso_simulation.h"
#include "gasphase/gasphasecv.h"
#include "particlephase/particlephase.h"
#include "particle.h"

#include <list>

class NetworkSimulation : public MesoSimulation {

protected:
	
	/// Simulation ID
	int ID;

	/// List of neightbourgs Simulations
	std::list<NetworkSimulation*> neighbourgs;

	/// List of buffers 
	std::list<std::list<std::shared_ptr<Particle>>> send_buffer;

	/// List of buffers 
	std::list<std::shared_ptr<Particle>> receive_buffer;

	/// List of concentrating species
	std::list<std::string> c_species;

	/// list of concentrating species molar fraction (ordered with respect of c_species)
	std::list<double> c_species_m_fraction;

	/// List of bath species
	std::list<std::string> b_species;

	//list of bath species molar fraction(ordered with respect of b_species)
	std::list<double> b_species_m_fraction;

	/// Pressure [Pa]
	double pressure;

	/// Start Temperature [K]
	double Temperature;

public:

	/// Blank Constructor
	NetworkSimulation();

	/// Construcor for Simulation Network
	///	\param: int _ID: Instance ID 
	///	\param: double _T: starting temperature [K]
	///	\param: double _p: pressure [pa]
	///	\param: std::list<std::string> c_species: Condensing species IUPAC names
	///	\param:	std::list<double> c_species_m_fraction: Condensing species molar fraction
	///	\param: std::list<std::string> b_species: Bath species IUPAC names
	///	\param:	std::list<double> b_species_m_fraction: Bath species molar fraction

	NetworkSimulation(int _ID, double _T, double _p,
		std::list<std::string> _c_species, std::list<double> _c_species_m_fraction,
		std::list<std::string> _b_species, std::list<double> _b_species_m_fraction);

	///	Criterion for selecting the particles to exchange with the neighbourgs
	virtual void select_particles() {}

	/// Ads a neighbourg to the simulation
	void add_neighbourg(NetworkSimulation* _n) { 
		neighbourgs.push_back(_n); 
		std::list<std::shared_ptr<Particle>> s_buff;
		send_buffer.push_back(s_buff);
	};

	/// Sends particles to neighbourgs simulations
	void send_particles();

	/// Gets particles from other simulation
	void get_particles(std::list<std::shared_ptr<Particle>> _part);

	/// Mesoscopic simulation 
	virtual void run_simulation(double _s_time, double _end_time) {};


};




#endif
