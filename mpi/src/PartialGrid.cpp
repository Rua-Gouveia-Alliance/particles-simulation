#include "PartialGrid.hpp"
#include "Mass.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <mpi.h>
#include <omp.h>

std::vector<long> PartialGrid::get_adjacent_cells(long ci, long ncside) {
  std::vector<long> adjacent;
  long cx = ci % ncside;
  long cy = ci / ncside;

  for (long ay = -1; ay <= 1; ay++) {
    for (long ax = -1; ax <= 1; ax++) {
      if (ay == 0 && ax == 0)
        continue;

      long nx = (cx + ax + ncside) % ncside;
      long ny = (cy + ay + ncside) % ncside;

      long idx = nx + ny * ncside;
      adjacent.push_back(idx);
    }
  }

  return adjacent;
}

std::vector<Mass> PartialGrid::_get_adjacent_masses(Cell &cell) {
  double side = cell.side;
  const std::vector<long> &adjacent_idx =
      get_adjacent_cells(cell.id(), _ncside);
  std::vector<Mass> adjacent_masses;

  for (long i : adjacent_idx) {
    PartialCell &adjacent_cell = adjacent_cells[i];
    Mass &mass = adjacent_cells[i].mass;

    // wrapping in x direction
    if (adjacent_cell.x >= cell.x + side * 2) {
      mass.x = mass.x - side;
    } else if (adjacent_cell.x < cell.x - side) {
      mass.x = mass.x + side;
    }

    // wrapping in y direction
    if (adjacent_cell.y >= cell.y + side * 2) {
      mass.y = mass.y - side;
    } else if (adjacent_cell.y < cell.y - side) {
      mass.y = mass.y + side;
    }

    adjacent_masses.push_back(mass);
  }

  return adjacent_masses;
}

long PartialGrid::_get_particle_index(Particle &p) {
  double temp_px, temp_py;

  // wrap particle around
  do {
    temp_px = p.x;
    temp_py = p.y;

    if (p.x < 0)
      p.x += _side;
    else if (p.x >= _side)
      p.x -= _side;

    if (p.y < 0)
      p.y += _side;
    else if (p.y >= _side)
      p.y -= _side;
  } while (temp_px != p.x || temp_py != p.y);

  // calculate cell index
  double cell_size = _side / _ncside;
  long cell_x = static_cast<long>(p.x / cell_size) % _ncside;
  long cell_y = static_cast<long>(p.y / cell_size) % _ncside;

  return cell_x + cell_y * _ncside;
}

void PartialGrid::_add_particle_to_cell(Particle &p) {
  long cell_idx = _get_particle_index(p);
  if (local_cells.find(cell_idx) != local_cells.end()) {
    local_cells[cell_idx].add_particle(p);
  } else {
    adjacent_cells[cell_idx].add_particle(p);
  }
}

void PartialGrid::_update_local_masses() {
  for (auto &it : local_cells) {
    it.second.update_mass();
  }
}

void PartialGrid::_update_local_cells() {
  std::vector<Particle> new_particles;
  long cell_count = local_cells.size();

  for (auto &it : local_cells) {
    Cell &cell = it.second;
    std::vector<Mass> masses = _get_adjacent_masses(cell);
    std::vector<Particle> part = cell.update_particles(adjacent_cells);
    new_particles.reserve(new_particles.size() + part.size());
    std::move(part.begin(), part.end(), std::back_inserter(new_particles));
  }

  _first_particle = _default_first_particle;
  for (auto &p : new_particles) {
    if (p.first_particle)
      _first_particle = p;
    _add_particle_to_cell(p);
  }

  for (auto &it : local_cells) {
    it.second.finish_update();
  }
}

void PartialGrid::update() {

  // for (const auto &id : _adjacent_ranks) {
  //   Mass mass = MPI_RECV_I(id);
  //   _adjacent_cells[mass.id()] = PartialCell(id, mass);
  // }

  _update_local_masses();

  // for (const auto &it : _local_cells) {
  //   for (const auto& id : it.second.adjacent_ranks) {
  //       MPI_SEND_I(it.second.mass, id);
  //   }
  // }

  // WAIT(_adjacent_cells);

  _update_local_cells();

  // for (const auto &it : _adjacent_cells) {
  //   MPI_SEND_I(it.second.owner(), it.second.particles);
  //   it.second.clear();
  // }

  // for (const auto &id : _adjacent_ranks) {
  //   std::vector<Particle> new_particles = MPI_RECV(id);
  //   for (auto &p : new_particles) {
  //     _add_particle_to_cell(p);
  //   }
  // }
  // MPI_WAIT_ALL(_adjacent_ranks);
}

void PartialGrid::sync_first_particle() {
  if (_rank == 0) {
    if (&_first_particle != &_default_first_particle)
      return;
    // _first_particle = MPI_RECV(ANY_RANK);
  } else if (&_first_particle != &_default_first_particle) {
    // MPI_SEND(0);
  }
}

void PartialGrid::print_cells() const {
  for (const auto &it : local_cells) {
    std::cout << "Cell " << it.second.id() << std::endl;
    it.second.print_particles();
    std::cout << std::endl;
  }
}

Particle PartialGrid::get_first_particle() const { return _first_particle; }

long PartialGrid::get_collisions() const {
  long total = 0;
  for (const auto &it : local_cells) {
    total += it.second.collisions;
  }
  return total;
}
