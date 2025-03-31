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

#define MASS_UPDATE 1
#define PARTICLE_UPDATE 2
#define FINAL_STATE 3

PartialGrid::PartialGrid(int rank, int max_rank, double side, long ncside)
    : _rank(rank), _max_rank(max_rank), _side(side), _ncside(ncside),
      _first_particle(_default_first_particle) {

  int m_count = 4;
  int m_blocklengths[4] = {1, 1, 1, 1};
  MPI_Datatype m_types[4] = {MPI_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint m_displacements[4] = {offsetof(Mass, id), offsetof(Mass, x),
                                 offsetof(Mass, y), offsetof(Mass, val)};
  MPI_Type_create_struct(m_count, m_blocklengths, m_displacements, m_types,
                         &mpi_mass_t);
  MPI_Type_commit(&mpi_mass_t);

  int p_count = 7;
  int p_blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
  MPI_Datatype p_types[7] = {MPI_INT,    MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL};
  MPI_Aint p_displacements[7] = {offsetof(Particle, id),
                                 offsetof(Particle, x),
                                 offsetof(Particle, y),
                                 offsetof(Particle, vx),
                                 offsetof(Particle, vy),
                                 offsetof(Particle, m),
                                 offsetof(Particle, first_particle)};
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
    Mass mass;
    double adjacent_cell_x, adjacent_cell_y;

    if (local_cells.find(i) != local_cells.end()) {
      Cell &adjacent_cell = local_cells.at(i);
      adjacent_cell_x = adjacent_cell.x;
      adjacent_cell_y = adjacent_cell.y;
      mass = adjacent_cell.mass;
    } else {
      PartialCell &adjacent_cell = adjacent_cells.at(i);
      adjacent_cell_x = adjacent_cell.x;
      adjacent_cell_y = adjacent_cell.y;
      mass = adjacent_cell.mass;
    }

    if (_rank == 2 && cell.id() == 8) {
      for (const auto &mass : adjacent_masses) {
        fprintf(stdout, "Mass %li ", mass.id);
        fprintf(stdout, "x: %#.6f ", mass.x);
        fprintf(stdout, "y: %#.6f ", mass.y);
        fprintf(stdout, "m: %#.6f\n", mass.val);
      }
    }

    // wrapping in x direction
    if (adjacent_cell_x >= cell.x + side * 2) {
      mass.x -= side;
    } else if (adjacent_cell_x < cell.x - side) {
      mass.x += side;
    }

    // wrapping in y direction
    if (adjacent_cell_y >= cell.y + side * 2) {
      mass.y -= side;
    } else if (adjacent_cell_y < cell.y - side) {
      mass.y += side;
    }

    adjacent_masses.push_back(mass);
  }

  if (_rank == 2 && cell.id() == 8) {
    for (const auto &mass : adjacent_masses) {
      fprintf(stdout, "Mass %li ", mass.id);
      fprintf(stdout, "x: %#.6f ", mass.x);
      fprintf(stdout, "y: %#.6f ", mass.y);
      fprintf(stdout, "m: %#.6f\n", mass.val);
    }
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
    std::vector<Particle> part = cell.update_particles(masses);
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

void PartialGrid::_print_trace() {
  if (_rank != 2)
    return;

  for (const auto &it : local_cells) {
    for (const auto &p : it.second.particles()) {
      fprintf(stdout, "Particle %lli: ", p.id);
      fprintf(stdout, "mass=%#.6f ", p.m);
      fprintf(stdout, "x=%#.6f ", p.x);
      fprintf(stdout, "y=%#.6f ", p.y);
      fprintf(stdout, "vx=%#.6f ", p.vx);
      fprintf(stdout, "vy=%#.6f\n", p.vy);
    }
  }

  for (long long i = 0; i < 9; ++i) {
    for (const auto &it : local_cells) {
      if (it.first == i) {
        fprintf(stdout, "Cell %lli ", i);
        fprintf(stdout, "x: %#.6f ", it.second.mass.x);
        fprintf(stdout, "y: %#.6f ", it.second.mass.y);
        fprintf(stdout, "m: %#.6f\n", it.second.mass.val);
        break;
      }
    }

    for (const auto &it : adjacent_cells) {
      if (it.first == i) {
        fprintf(stdout, "Cell %lli ", i);
        fprintf(stdout, "x: %#.6f ", it.second.mass.x);
        fprintf(stdout, "y: %#.6f ", it.second.mass.y);
        fprintf(stdout, "m: %#.6f\n", it.second.mass.val);
        break;
      }
    }
  }
}

void PartialGrid::update() {
  long size = _adjacent_ranks.size();
  std::vector<MPI_Request> requests(size * 2);
  std::vector<std::vector<Mass>> masses;
  masses.reserve(size);

  long i = 0;
  // Receiving mass updates from each adjacent rank
  for (const auto &it : _adjacent_ranks) {
    int rank = it.first;
    long element_count = it.second;
    masses.push_back(std::vector<Mass>(element_count));
    MPI_Irecv(masses.back().data(), element_count, mpi_mass_t, rank,
              MASS_UPDATE, MPI_COMM_WORLD, &requests[i++]);
  }

  _update_local_masses();

  // Calculating the mass updates to send to each adjacent rank
  std::unordered_map<int, std::vector<Mass>> mass_updates;
  mass_updates.reserve(size);
  for (const auto &it : local_cells) {
    for (const auto &id : it.second.adjacent_ranks) {
      mass_updates[id].push_back(it.second.mass);
    }
  }

  // Sending the mass updates to each adjacent rank
  for (auto &it : mass_updates) {
    MPI_Isend(it.second.data(), it.second.size(), mpi_mass_t, it.first,
              MASS_UPDATE, MPI_COMM_WORLD, &requests[i++]);
  }

  std::vector<MPI_Status> statuses(size * 2);
  // Waiting for mass updates isend/irecv
  MPI_Waitall(size * 2, requests.data(), statuses.data());
  for (const auto &array : masses) {
    for (const auto &mass : array) {
      adjacent_cells.at(mass.id).mass = mass;
    }
  }

  _print_trace();

  _update_local_cells();

  // Calculating the particle updates to send to each adjacent rank
  std::unordered_map<int, std::vector<Particle>> particle_updates;
  particle_updates.reserve(size);
  for (auto &it : adjacent_cells) {
    PartialCell &cell = it.second;
    std::vector<Particle> &update = particle_updates[cell.owner];
    std::move(cell.particles.begin(), cell.particles.end(),
              std::back_inserter(update));
    cell.clear();
  }

  i = 0;
  // Sending the particle updates to each adjacent rank
  for (auto &it : particle_updates) {
    MPI_Isend(it.second.data(), it.second.size(), mpi_particle_t, it.first,
              PARTICLE_UPDATE, MPI_COMM_WORLD, &requests[i++]);
  }

  // Receiving particle updates from each adjacent rank
  for (const auto &it : _adjacent_ranks) {
    int count;
    MPI_Status status;
    MPI_Probe(it.first, PARTICLE_UPDATE, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, mpi_particle_t, &count);

    std::vector<Particle> particles(count);
    MPI_Recv(particles.data(), count, mpi_particle_t, it.first, PARTICLE_UPDATE,
             MPI_COMM_WORLD, &status);

    for (auto &p : particles) {
      long idx = _get_particle_index(p);
      local_cells.at(idx).add_updated_particle(p);
    }
  }

  for (auto &it : local_cells) {
    it.second.check_collisions();
  }

  // Waiting for particle updates isend
  MPI_Waitall(size, requests.data(), statuses.data());
}

void PartialGrid::sync_final_state() {
  typedef struct {
    Particle particle;
    long collisions;
  } final_state_t;

  MPI_Datatype mpi_final_state_t;
  int count = 2;
  int blocklengths[2] = {1, 1};
  MPI_Datatype types[2] = {mpi_particle_t, MPI_LONG};
  MPI_Aint displacements[2] = {offsetof(final_state_t, particle),
                               offsetof(final_state_t, collisions)};
  MPI_Type_create_struct(count, blocklengths, displacements, types,
                         &mpi_final_state_t);
  MPI_Type_commit(&mpi_final_state_t);

  if (_rank == 0) {
    MPI_Status status;
    final_state_t state;

    for (int i = 1; i < _max_rank + 1; ++i) {
      MPI_Recv(&state, 1, mpi_final_state_t, i, FINAL_STATE, MPI_COMM_WORLD,
               &status);
      _remote_collisions += state.collisions;
      if (state.particle.first_particle)
        _first_particle = state.particle;
    }
  } else {
    final_state_t state = {_first_particle, get_collisions()};
    MPI_Send(&state, 1, mpi_final_state_t, 0, FINAL_STATE, MPI_COMM_WORLD);
  }

  MPI_Type_free(&mpi_mass_t);
  MPI_Type_free(&mpi_particle_t);
  MPI_Type_free(&mpi_final_state_t);
}

Particle PartialGrid::get_first_particle() const { return _first_particle; }

long PartialGrid::get_collisions() const {
  long total = 0;
  for (const auto &it : local_cells) {
    total += it.second.collisions;
  }
  return total + _remote_collisions;
}
