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

#ifndef CONSTRAINEDLANGEVINPARTICLEPHASE_H
#define CONSTRAINEDLANGEVINPARTICLEPHASE_H

#include "dynamicparticlephase.h"
#include "pbmfractalparticlephase.h"

#include <time.h>

//#define OMP
#include <omp.h>

template<typename A>
class ConstrainedLangevinParticlePhase : public DynamicParticlePhase<A> {

	/// Function for calculating the value of a 3d index a grid in a periodic domain
	///	\param	int c_idx: cell actual component position
	///	\param int SIDE: Number of cells in the dimension
	///	\param double& shift: shift to impose in the point coordinates at the correspondig component
	///	\param double domain_side: lenght of the control volume side
	int periodic_position(int c_idx, int SIDE, double& shift, double domain_side);

public:
	
	/// Constructor
    ConstrainedLangevinParticlePhase(double _volume):
        DynamicParticlePhase<A>(_volume){}

	/// Constructor from PBM Particle phase
	///	\param PBMFractalParticlePhase<PBMAggregate<Particle>>& _fractal_pp PBM particle phase
	ConstrainedLangevinParticlePhase(PBMFractalParticlePhase<PBMAggregate<Particle>>& _fractal_pp, double _volume, double T) :
		DynamicParticlePhase<A>(_volume) 
	{

		double PBM_Volume = _fractal_pp.get_volume();

		int N = _fractal_pp.get_aggregates_number();

		//if volumes are the same

		if (PBM_Volume == _volume) {
			// Spawn aggregates from PBM aggregates
			for (int i = 0; i < N; i++) {
				std::shared_ptr<PBMAggregate<Particle>> pbm_agg = _fractal_pp.get_aggregate(i);
				this->add_particle_lc(pbm_agg->get_n_monomers(), pbm_agg->get_species());

			}
		}
		else if (_volume < PBM_Volume) {

			// Choose ramdom aggregates from the PBM Particle Phase
			double pbm_p_density = N / PBM_Volume;
			int n_particles = (N*_volume) / PBM_Volume;
			std::cout << "PBM Particles:" << N << std::endl;
			std::cout << "Langevin Particles" << n_particles << std::endl;

			srand(time(NULL));
			int count = 0;
			while (count < n_particles) {
				int rand_index = rand() % N;
				std::shared_ptr<PBMAggregate<Particle>> pbm_agg = _fractal_pp.get_aggregate(rand_index);
				this->add_particle_lc(pbm_agg->get_n_monomers(), pbm_agg->get_species());
				count++;
			}

		}
		else { // if the volume is bigger than the PBM one
			// TO BE IMPLEMENTED
		}

		// Initialize velocities
		this->initialize_velocities(T);
	}

	/// Constructor for Particle Phase starting from a list of aggregates (hotstart/recovery)
	ConstrainedLangevinParticlePhase(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& _aggregates, double _volume, double _T):
		DynamicParticlePhase<A>(_volume) {

		// Assign the list
		this->aggregates = _aggregates;

		// Initialize velocities
		this->initialize_velocities(_T);
		
		// Displace aggregates in the grid
		this->grid.displace_aggregates(this->aggregates);
	}

    /// Saves VTK data
	///	\param int iteration: iteration in which the snapshot is saved
	///	\param s_idx: int streamline index
    void save_vtk(int iteration, std::string root_path);

    /// Time step
	///	\param double dt: timestep [sec]
	///	\param const GasPhase& gp: GasPhase object [by reference]
	///	\param const NucleationTheory& nt: Nucleation Theory object [by reference]
	///	\param double& T: Temperature[K] [by reference]
    double timestep(double dt, const GasPhase& gp, const NucleationTheory& nt, double& T);

	/// Linked-Cell time step
	///	\param double dt: timestep [sec]
	///	\param const GasPhase& gp: GasPhase object [by reference]
	///	\param const NucleationTheory& nt: Nucleation Theory object [by reference]
	///	\param double& T: Temperature[K] [by reference]
	double timestep_lc(double dt, const GasPhase& gp, const NucleationTheory& nt, double& T);

	/// Returns cardinality of the active particles in the simulation
    int get_aggregates_cardinality();

	/// Updates graph constraints values
	void update_constraints();

	/// biggest spherical enclosure diameter
	double get_biggest_spherical_enclosure(int& _id) const;


protected:
	
	/// Motion
	///	\param dt timestep [s]
	///	\param gp Gas Phase object
	///	\param T Temperature [K]
	void motion(double dt, const GasPhase& gp, double T);

	/// Motion (Linked Cell Version)
	///	\param dt timestep [s]
	///	\param gp Gas Phase object
	///	\param T Temperature [K]
	void motion_lc(double dt, const GasPhase& gp, double T);

    ///Coagulation
    void coagulation();

	/// Coagulation (Linked Cell Version)
	void coagulation_lc();

    /// Creates a Particle in the simulation
    void add_particle(double n, Species s);

	/// Creates a Particle in the simulation (Linked Cell Version)
	void add_particle_lc(double n, Species s);

    /// Checks if some particles has coalesced
    void rearrange_aggregate();

	/// Keeps constant the number of aggregates in the simulation
	void aggregates_number_balance();

	/// Initialize Velocity
	///	\param double T Temperature [K]
	void initialize_velocities(double T);

    /// Evaluate nucleation events for species s and return the molecules consumption
    /// [#/m3 s]
    /// \param dt timestep [s]
    /// \param J particle nucleation rate [#/m3 s]
    /// \param s condensing species
    double nucleation(double dt, double J, double j, Species s);

	/// Evaluate nucleation events for species s and return the molecules consumption (Linked Cell Version)
	/// [#/m3 s]
	/// \param dt timestep [s]
	/// \param J particle nucleation rate [#/m3 s]
	/// \param s condensing species
	double nucleation_lc(double dt, double J, double j, Species s);

	/// Expands the volume and shifts the aggregates
	/// Correct volume due to gas expansion
	/// \param dt timestep [s]
	/// \param gp gas phase surrounding particles
	void volume_resize(double dt, const GasPhase& gp);


};

#include "constrainedlangevinparticlephase.cpp"

#endif // CONSTRAINEDLANGEVINPARTICLEPHASE_H
