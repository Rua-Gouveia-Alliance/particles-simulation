#pragma once

#include "Cell.hpp"
#include "Mass.hpp"
#include "PartialCell.hpp"
#include "Particle.hpp"
#include <mpi.h>
#include <unordered_map>
#include <vector>

class PartialGrid {
private:
  int _rank;
  double _side;
  long _ncside;
  Particle &_first_particle;
  Particle _default_first_particle = Particle();
  std::vector<int> _adjacent_ranks;
  MPI_Datatype mpi_mass_t;
  MPI_Datatype mpi_particle_t;

  std::vector<Mass> _get_adjacent_masses(Cell &cell);
  long _get_particle_index(Particle &p);
  void _add_particle_to_cell(Particle &p);
  void _update_local_masses();
  void _update_local_cells();

public:
  std::unordered_map<int, PartialCell> adjacent_cells;
  std::unordered_map<int, Cell> local_cells;

  PartialGrid(int rank, double side, long ncside);

  static std::vector<long> get_adjacent_cells(long ci, long ncside);
  void add_local_cell(Cell cell) { local_cells.emplace(cell.id(), cell); }
  void add_adjacent_cell(PartialCell cell) {
    adjacent_cells.emplace(cell.id(), cell);
  }
  void add_adjacent_rank(long rank) { _adjacent_ranks.push_back(rank); }

  void update();
  void sync_first_particle();

  void print_cells() const;
  Particle get_first_particle() const;
  long get_collisions() const;
};
