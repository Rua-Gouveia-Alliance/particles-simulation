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
  double _side;
  int _rank, _max_rank;
  long _ncside, _remote_collisions = 0;
  Particle &_first_particle;
  Particle _default_first_particle = Particle();
  std::unordered_map<int, long> _adjacent_ranks;
  MPI_Datatype mpi_mass_t, mpi_particle_t;

  std::vector<Mass> _get_adjacent_masses(Cell &cell);
  long _get_particle_index(Particle &p);
  void _add_particle_to_cell(Particle &p);
  void _update_local_masses();
  void _update_local_cells();

public:
  std::unordered_map<int, PartialCell> adjacent_cells;
  std::unordered_map<int, Cell> local_cells;

  PartialGrid(int rank, int max_rank, double side, long ncside);

  static std::vector<long> get_adjacent_cells(long ci, long ncside);
  void add_local_cell(Cell cell) { local_cells.try_emplace(cell.id(), cell); }
  void add_adjacent_cell(PartialCell cell) {
    adjacent_cells.try_emplace(cell.id(), cell);
  }
  void add_adjacent_rank(int rank) { _adjacent_ranks.emplace(rank, 0); }
  void increment_adjacent_rank(int rank) { ++_adjacent_ranks.at(rank); }

  void update();
  void sync_final_state();

  void print_cells() const;
  Particle get_first_particle() const;
  long get_collisions() const;
};
