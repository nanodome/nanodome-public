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

#ifndef CONSTRAINER_H
#define CONSTRAINER_H

#include "tetra_vertex.h"
#include "tetra_constraint.h"
#include "tetra_bond.h"
#include "tetra_tetrahedron.h"
#include "tetra_link.h"
#include "tetra_constants.h"

#include "dynamicparticle.h"
#include "particlebond.h"


#include <vector>
#include <list>
#include <valarray>
#include <string>
#include <random>
#include <map>

//class DynamicParticle;

class Constrainer
{
    /// Mapping among local and global indices
    std::map<int, int> global_to_local;

    /// Counter for vertices local index
    int counter;

    /// Vector of the vertices composing the aggregate
    std::vector<Vertex> vertices;

    /// Adjacency matrix for the connections among vertices
    std::vector<std::vector<PhysicalLink > > bond_edges;

    /// Adjacency matrix for the contraints among vertices
    std::vector<std::vector<Constraint > > constraint_edges;

    /// List of tetrahedra composing the structure
    std::list<Tetrahedron> tetrahedra;

    /*----- Private Methods -----*/

    /// Creates contraints for a 3 vertices graph
    void constraints3();

    /// Creates constraints for a 4 vertices graph (tetrahedron)
    /// \param: std::vector<Vertex*> vertices list of vertices for creating a thetrahedron
    void constraints4(std::vector<Vertex*> vertices);

    /// Creates constraints for a graph with |V| > 4 vertices
    void constraintsN();

    /// Breadth First Search to find a connected component of the graph
    /// for creating the first tetrahedron
    /// \param index: index of the starting vertex for the search
    std::vector<Vertex*> BFS(int _index);


public:

    /// Blank Constructor
    Constrainer();

    /// Constructor
    /// \param DynamicParticle _p: first aggregate particle
    Constrainer(std::shared_ptr< DynamicParticle >_p);

    /// Create constraints
    void create_constraints();

    /// Returns number of vertices
    int get_n_vert() const { return vertices.size(); }

    /// Returns the vertices
    std::vector<Vertex> get_vertices()const { return vertices; }

    /// Returns a list with the edges
    std::list<Link*> get_edges();

	/// Returns a vector with couples of global indices for each edge (constraint and physical) and the type
	/// triple i, j, type
	std::vector<std::tuple<int, int, int>> get_edges_indices();

    /// Returns all the edges by pairs of estremes (local indices)
    std::vector<std::pair<int, int>> get_edge_extremes();

    /// Returns a list with the edges types
    std::list<std::string> get_edges_types();

    /// Returns number of edges
    int get_n_edges();

    /// Reset the Graph
    void reset_graph();

    /// Shake algorithm
    bool SHAKE();

    /// Adds Particles to the graph
    void add_particles(std::list<std::shared_ptr<DynamicParticle>> p_list,
                       std::list < std::shared_ptr<ParticleBond<DynamicParticle>>> b_list);

    /// Updates bond distances
    void update_bonds(std::list < std::shared_ptr<ParticleBond<DynamicParticle>>> b_list);

};
#endif //CONSTRAINER_H
