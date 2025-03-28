#pragma once

#include "Cell.hpp"
#include "Mass.hpp"
#include "PartialCell.hpp"
#include "Particle.hpp"
#include <unordered_map>
#include <vector>

class PartialGrid {
private:
  double _side;
  long _ncside;
  Particle &_first_particle;
  Particle _default_first_particle = Particle();
  std::vector<int> _adjacent_ranks;

  std::vector<Mass> _get_adjacent_masses(Cell &cell);
  long _get_particle_index(Particle &p);
  void _add_particle_to_cell(Particle &p);
  void _update_local_masses();
  void _update_local_cells();

public:
  std::unordered_map<int, PartialCell> adjacent_cells;
  std::unordered_map<int, Cell> local_cells;

  PartialGrid(double side, long ncside, std::vector<int> adjacent_ranks)
      : _side(side), _ncside(ncside), _adjacent_ranks(adjacent_ranks),
        _first_particle(_default_first_particle) {}

  static std::vector<long> get_adjacent_cells(long ci, long ncside);
  void add_local_cell(Cell &cell) { local_cells[cell.id()] = cell; }
  void add_adjacent_cell(PartialCell &cell) {
    adjacent_cells[cell.id()] = cell;
  }

  void update();

  void print_cells() const;
  Particle get_first_particle() const;
  long get_collisions() const;
};
