#include "PartialGrid.hpp"
#include "Mass.hpp"
#include "Particle.hpp"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <mpi.h>
#include <omp.h>
#include <stddef.h>
#include <unordered_map>
#include <vector>

#define FIRST_PARTICLE 1
#define MASS_UPDATE 2

PartialGrid::PartialGrid(int rank, double side, long ncside)
    : _rank(rank), _side(side), _ncside(ncside),
      _first_particle(_default_first_particle) {

  int m_count = 4;
  int m_blocklengths[4] = {1, 1, 1, 1};
  MPI_Datatype m_types[4] = {MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint m_displacements[4] = {offsetof(Mass, id), offsetof(Mass, x),
                                 offsetof(Mass, y), offsetof(Mass, val)};
  MPI_Type_create_struct(m_count, m_blocklengths, m_displacements, m_types,
                         &mpi_mass_t);
  MPI_Type_commit(&mpi_mass_t);

  int p_count = 6;
  int p_blocklengths[6] = {1, 1, 1, 1, 1, 1};
  MPI_Datatype p_types[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL};
  MPI_Aint p_displacements[6] = {
      offsetof(Particle, x),  offsetof(Particle, y),
      offsetof(Particle, vx), offsetof(Particle, vy),
      offsetof(Particle, m),  offsetof(Particle, first_particle)};
  MPI_Type_create_struct(p_count, p_blocklengths, p_displacements, p_types,
                         &mpi_particle_t);
  MPI_Type_commit(&mpi_particle_t);
}

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
    PartialCell &adjacent_cell = adjacent_cells.at(i);
    Mass &mass = adjacent_cell.mass;

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
    local_cells.at(cell_idx).add_particle(p);
  } else {
    adjacent_cells.at(cell_idx).add_particle(p);
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
  long size = _adjacent_ranks.size();
  std::vector<MPI_Request> requests(size * 2);
  std::vector<std::vector<Mass>> masses;
  masses.reserve(size);

  long i = 0;
  for (const auto &it : _adjacent_ranks) {
    int rank = it.first;
    long element_count = it.second;
    masses.push_back(std::vector<Mass>(element_count));
    MPI_Irecv(masses.back().data(), element_count, mpi_mass_t, rank,
              MASS_UPDATE, MPI_COMM_WORLD, &requests[i++]);
  }

  _update_local_masses();

  std::unordered_map<int, std::vector<Mass>> updates;
  updates.reserve(size);
  for (const auto &it : local_cells) {
    for (const auto &id : it.second.adjacent_ranks) {
      updates[id].push_back(it.second.mass);
    }
  }

  for (auto &it : updates) {
    MPI_Isend(it.second.data(), it.second.size(), mpi_mass_t, it.first,
              MASS_UPDATE, MPI_COMM_WORLD, &requests[i++]);
  }

  std::vector<MPI_Status> statuses(size * 2);
  MPI_Waitall(size * 2, requests.data(), statuses.data());
  for (const auto &array : masses) {
    for (const auto &mass : array) {
      adjacent_cells.at(mass.id).mass = mass;
    }
  }

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

    MPI_Status status;
    MPI_Recv(&_first_particle, 1, mpi_particle_t, MPI_ANY_SOURCE,
             FIRST_PARTICLE, MPI_COMM_WORLD, &status);
  } else if (&_first_particle != &_default_first_particle) {
    MPI_Send(&_first_particle, 1, mpi_particle_t, 0, FIRST_PARTICLE,
             MPI_COMM_WORLD);
  }

  MPI_Type_free(&mpi_mass_t);
  MPI_Type_free(&mpi_particle_t);
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
