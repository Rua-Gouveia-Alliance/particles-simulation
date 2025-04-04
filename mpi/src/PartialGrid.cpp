#include "PartialGrid.hpp"
#include "Mass.hpp"
#include "Particle.hpp"
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <unordered_map>
#include <vector>

#define MASS_UPDATE 1
#define PARTICLE_UPDATE 2
#define FINAL_STATE 3

PartialGrid::PartialGrid(int rank, int nprocs, double side, long ncside)
    : _rank(rank), _nprocs(nprocs), _side(side), _ncside(ncside),
      _first_particle({.first_particle = false}) {

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

void PartialGrid::add_local_cell(Cell &cell) {
  if (cell.adjacent_ranks.size() == 0) {
    fully_local_cells.try_emplace(cell.id(), cell);
  } else {
    partially_local_cells.try_emplace(cell.id(), cell);
  }
}

std::vector<Mass> PartialGrid::_get_adjacent_masses(Cell &cell) {
  double side = cell.side;
  const std::vector<long> &adjacent_idx =
      get_adjacent_cells(cell.id(), _ncside);
  std::vector<Mass> adjacent_masses;

  for (long i : adjacent_idx) {
    Mass mass;
    double adjacent_cell_x, adjacent_cell_y;

    if (fully_local_cells.find(i) != fully_local_cells.end()) {
      Cell &adjacent_cell = fully_local_cells.at(i);
      adjacent_cell_x = adjacent_cell.x;
      adjacent_cell_y = adjacent_cell.y;
      mass = adjacent_cell.mass;
    } else if (partially_local_cells.find(i) != partially_local_cells.end()) {
      Cell &adjacent_cell = partially_local_cells.at(i);
      adjacent_cell_x = adjacent_cell.x;
      adjacent_cell_y = adjacent_cell.y;
      mass = adjacent_cell.mass;
    } else {
      PartialCell &adjacent_cell = adjacent_cells.at(i);
      adjacent_cell_x = adjacent_cell.x;
      adjacent_cell_y = adjacent_cell.y;
      mass = adjacent_cell.mass;
    }

    // wrapping in x direction
    if (adjacent_cell_x >= cell.x + side * 2) {
      mass.x = mass.x - _side;
    } else if (adjacent_cell_x < cell.x - side) {
      mass.x = mass.x + _side;
    }

    // wrapping in y direction
    if (adjacent_cell_y >= cell.y + side * 2) {
      mass.y = mass.y - _side;
    } else if (adjacent_cell_y < cell.y - side) {
      mass.y = mass.y + _side;
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
  long idx = _get_particle_index(p);
  if (fully_local_cells.find(idx) != fully_local_cells.end()) {
    fully_local_cells.at(idx).add_particle(p);
  } else if (partially_local_cells.find(idx) != partially_local_cells.end()) {
    partially_local_cells.at(idx).add_particle(p);
  } else {
    adjacent_cells.at(idx).add_particle(p);
  }
}

void PartialGrid::_update_local_masses(std::unordered_map<int, Cell> &cells) {
  for (auto &it : cells) {
    it.second.update_mass();
  }
}

void PartialGrid::_update_local_cells(std::unordered_map<int, Cell> &cells) {
  std::vector<Particle> new_particles;
  long cell_count = cells.size();

  const auto &begin = cells.begin();
#pragma omp for schedule(dynamic)
  for (long i = 0; i < cell_count; ++i) {
    Cell &cell = std::next(begin, i)->second;
    std::vector<Mass> masses = _get_adjacent_masses(cell);
    std::vector<Particle> part = cell.update_particles(masses);

#pragma omp critical
    {
      new_particles.reserve(new_particles.size() + part.size());
      std::move(part.begin(), part.end(), std::back_inserter(new_particles));
    }
  }

  for (auto &p : new_particles)
    _add_particle_to_cell(p);
}

void PartialGrid::_finish_local_updates() {
  const auto &begin_0 = fully_local_cells.begin();
  long fl_size = fully_local_cells.size();
#pragma omp for
  for (long i = 0; i < fl_size; i++) {
    const auto &it = std::next(begin_0, i);
    it->second.finish_update();
  }

  const auto &begin_1 = partially_local_cells.begin();
  long pl_size = partially_local_cells.size();
#pragma omp for
  for (long i = 0; i < pl_size; i++) {
    const auto &it = std::next(begin_1, i);
    it->second.finish_update();
  }
}

void PartialGrid::update(long long time_steps) {
  // This rank is not responsible for any cells
  if (fully_local_cells.size() == 0 && partially_local_cells.size() == 0) {
    return;
  }

  long size = _adjacent_ranks.size();
  std::vector<MPI_Request> requests(size * 2);
  std::vector<std::vector<Mass>> masses;
  std::unordered_map<int, std::vector<Particle>> particle_updates;
  std::unordered_map<int, std::vector<Mass>> mass_updates;

  masses.reserve(size);
  particle_updates.reserve(size);
  mass_updates.reserve(size);

#pragma omp parallel
  for (long long ll = 0; ll < time_steps; ll++) {
    long i = 0;

    // Receiving mass updates from each adjacent rank
    for (const auto &it : _adjacent_ranks) {
      int rank = it.first;
      long element_count = it.second;
      masses.push_back(std::vector<Mass>(element_count));
      MPI_Irecv(masses.back().data(), element_count, mpi_mass_t, rank,
                MASS_UPDATE, MPI_COMM_WORLD, &requests[i++]);
    }

    _update_local_masses(partially_local_cells);

    // Calculating the mass updates to send to each adjacent rank
    mass_updates.clear();
    for (const auto &it : partially_local_cells) {
      for (const auto &id : it.second.adjacent_ranks) {
        mass_updates[id].push_back(it.second.mass);
      }
    }

    // Sending the mass updates to each adjacent rank
    const auto &begin_0 = mass_updates.begin();
#pragma omp for
    for (long j = 0; j < size; ++j) {
      const auto &it = std::next(begin_0, j);
      MPI_Isend(it->second.data(), it->second.size(), mpi_mass_t, it->first,
                MASS_UPDATE, MPI_COMM_WORLD, &requests[i++]);
    }

    _update_local_masses(fully_local_cells);

    std::vector<MPI_Status> statuses(size * 2);
    // Waiting for mass updates isend/irecv
    MPI_Waitall(size * 2, requests.data(), statuses.data());
    for (const auto &array : masses) {
#pragma omp for
      for (const auto &mass : array) {
        adjacent_cells.at(mass.id).mass = mass;
      }
    }

    _update_local_cells(partially_local_cells);

    // Calculating the particle updates to send to each adjacent rank
    particle_updates.clear();
    for (auto &it : adjacent_cells) {
      PartialCell &cell = it.second;
      std::vector<Particle> &update = particle_updates[cell.owner];
      std::move(cell.particles.begin(), cell.particles.end(),
                std::back_inserter(update));
      cell.clear();
    }

    // Sending the particle updates to each adjacent rank
    i = 0;
    const auto &begin_1 = particle_updates.begin();
#pragma omp for
    for (long j = 0; j < size; ++j) {
      const auto &it = std::next(begin_1, j);
      MPI_Isend(it->second.data(), it->second.size(), mpi_particle_t, it->first,
                PARTICLE_UPDATE, MPI_COMM_WORLD, &requests[i++]);
    }

    _update_local_cells(fully_local_cells);
    _finish_local_updates();

    const auto &begin_2 = fully_local_cells.begin();
    long fl_size = fully_local_cells.size();
#pragma omp for schedule(dynamic)
    for (long j = 0; j < fl_size; ++j) {
      const auto &it = std::next(begin_2, j);
      it->second.check_collisions();
    }

    // Receiving particle updates from each adjacent rank
    for (const auto &it : _adjacent_ranks) {
      int count;
      MPI_Status status;
      MPI_Probe(it.first, PARTICLE_UPDATE, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, mpi_particle_t, &count);

      std::vector<Particle> particles(count);
      MPI_Recv(particles.data(), count, mpi_particle_t, it.first,
               PARTICLE_UPDATE, MPI_COMM_WORLD, &status);

      for (auto &p : particles) {
        long idx = _get_particle_index(p);
        partially_local_cells.at(idx).add_updated_particle(p);
      }
    }

    const auto &begin_3 = partially_local_cells.begin();
    long pl_size = partially_local_cells.size();
#pragma omp for schedule(dynamic)
    for (long j = 0; j < pl_size; j++) {
      const auto &it = std::next(begin_3, j);
      it->second.check_collisions();
    }

    // Waiting for particle updates isend
    MPI_Waitall(size, requests.data(), statuses.data());
  }
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

  for (const auto &it : partially_local_cells) {
    for (const auto &p : it.second.particles) {
      if (p.first_particle)
        _first_particle = p;
    }
  }

  for (const auto &it : fully_local_cells) {
    for (const auto &p : it.second.particles) {
      if (p.first_particle)
        _first_particle = p;
    }
  }

  final_state_t state = {_first_particle, get_collisions()};
  std::vector<final_state_t> states(_nprocs);
  MPI_Gather(&state, 1, mpi_final_state_t, states.data(), 1, mpi_final_state_t,
             0, MPI_COMM_WORLD);
  if (_rank == 0) {
    for (const auto &s : states) {
      _remote_collisions += s.collisions;
      if (s.particle.first_particle) {
        _first_particle = s.particle;
      }
    }
  }

  MPI_Type_free(&mpi_mass_t);
  MPI_Type_free(&mpi_particle_t);
  MPI_Type_free(&mpi_final_state_t);
}

long PartialGrid::get_collisions() const {
  long total = 0;

  for (const auto &it : fully_local_cells) {
    total += it.second.collisions;
  }
  for (const auto &it : partially_local_cells) {
    total += it.second.collisions;
  }

  return total + _remote_collisions;
}
