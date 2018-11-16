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

#ifndef GRID_H
#define GRID_H

#include "cell.h"

//#include "aggregate/spatialaggregate.h"
#include "aggregate/rattleaggregate.h"
#include "dynamicparticle.h"
#include "i_o.h"

#include <vector>


class Grid {

	
	/// Linked-Cell grid
	std::vector<Cell<RATTLEAggregate<DynamicParticle>>> grid;

	/// Cell side (cubic)
	double cell_side;

	/// Number of cells for side
	int C;

	/// Control Volume side
	double volume_side;

	/// Grid side
	double grid_side;

	/// Number of cells
	int n_cells;

	/// position shift from centre of the simulation(centre of the control volume) to centre of the grid (leftmost vertex)
	double shift;

	///	Create the grid
	void grid_init();

	/// Erase the grid
	void grid_erase();

	/// Returns linear index in the grid (x, y, z)
	int get_linear(const int x1, const int x2, const int x3) {
		return x3*C*C + x2*C + x1;}

public:

	/// Blank Constructor
	Grid();

	/// Constructor
	///	\param double _volume Control Volume dimension [m3]
	/// \param int _dimension starting number of cell for side(eg. _dimension = 5 -> grid of 5x5x5 cells)
	Grid(double _volume, int _dimension);

	/// Grid Update/Creation
	///	\param std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble: particles ensamble
	///	\param double _V: Control Volume dimension [m3]
	///	\param double _agg_max_diameter: Max aggregate diameter [m]
	///	\param int _max_id: id of the biggest particle in the ensamble [index]
	void grid_update(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble,
					 double _V, double _agg_max_diameter, int _max_id);

	/// Aggregates Dispacement on the grid
	///	\param std::list<std::shared_ptr<SpatialAggregate<DynamicParticle>>> _ensamble: list of aggregates in the simulation
	void displace_aggregates(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble);

	/// Adds aggregates to the grid
	///	\param std::shared_ptr<SpatialAggregate<DynamicParticle>> _p0: pointer to the aggregate to add
	void insert_aggregate(std::shared_ptr<RATTLEAggregate<DynamicParticle>> _p0 );

	/// Moves aggregates on the grid
	///	\param std::shared_ptr<RATTLEAggregate<DynamicParticle>> _p0 pointer to the aggregate to move
	void move_aggregates();

	/// Returns total grid cells number
	int get_n_cells() const { return n_cells; }

	/// Returns grid number of cells for side
	int get_C() const { return C; }

	/// Returns the cell side of a single cell in the grid
	double get_cell_side() const { return cell_side; }

	/// Returns the grid side
	double get_grid_side()const { return grid_side; }

	/// Returns the cell of the relative coordinates
	///	\param std::valarray<double> _x point coordinates to retrieve the cell index (_x must be a valid set of coordinates in the domain)
	std::valarray<int> get_cell_index(std::valarray<double> _x) const;

	/// Overload operator "()" for access the grid with multiple indexes (3D)
	Cell<RATTLEAggregate<DynamicParticle>> &operator()(const int _x1, const int _x2, const int _x3); //(x, y, z)

	/// Overload operator "[]" for access the grid with single index
	Cell<RATTLEAggregate<DynamicParticle>> &operator[](const int lin_index);

};

#endif

