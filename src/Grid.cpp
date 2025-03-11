#include "Grid.hpp"
#include "Particle.hpp"
#include <iostream>

#include <omp.h>

Grid::Grid(double side, long ncside)
    : _side(side), _ncside(ncside), _first_particle(_default_first_particle) {}

void Grid::_add_particle_to_cell(Particle &p) {
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

  long cell_idx = cell_x + cell_y * _ncside;
#pragma omp critical
  _cells[cell_idx].add_particle(p);
}

void Grid::add_cell(Cell &c) { _cells.push_back(c); }

std::vector<long> Grid::get_adjacent_cells(long ci) {
  std::vector<long> adjacent;
  long cx = ci % _ncside;
  long cy = ci / _ncside;

  for (long ay = -1; ay <= 1; ay++) {
    for (long ax = -1; ax <= 1; ax++) {
      if (ay == 0 && ax == 0)
        continue;

      long nx = (cx + ax + _ncside) % _ncside;
      long ny = (cy + ay + _ncside) % _ncside;

      long idx = nx + ny * _ncside;
      adjacent.push_back(idx);
    }
  }

  return adjacent;
}

void Grid::update_cells() {
  long num_cells = _cells.size();
  std::vector<std::vector<Particle>> new_particles(num_cells);

#pragma omp parallel for
  for (long i = 0; i < num_cells; i++) {
    Cell &curr_cell = _cells[i];
    double cell_side;

    const std::vector<long> &adjacent_idx = get_adjacent_cells(i);
    std::vector<Cell> adjacent_cells;
    for (long j : adjacent_idx) {
      Cell new_cell = _cells[j];
      cell_side = curr_cell.side;

      // wrapping in x direction
      if (new_cell.x >= curr_cell.x + cell_side * 2) {
        new_cell.center_of_mass_x -= _side;
      } else if (new_cell.x < curr_cell.x - cell_side) {
        new_cell.center_of_mass_x += _side;
      }

      // wrapping in y direction
      if (new_cell.y >= curr_cell.y + cell_side * 2) {
        new_cell.center_of_mass_y -= _side;
      } else if (new_cell.y < curr_cell.y - cell_side) {
        new_cell.center_of_mass_y += _side;
      }

      adjacent_cells.push_back(new_cell);
    }

    new_particles[i] = curr_cell.update_particles(adjacent_cells);
  }

#pragma omp parallel for
  for (auto &v : new_particles) {
    for (auto &p : v) {
      if (p.first_particle)
        _first_particle = p;

      _add_particle_to_cell(p);
    }
  }

#pragma omp parallel for
  for (auto &c : _cells) {
    c.finish_update();
  }
}

void Grid::print_cells() const {
  for (long i = 0; i < _cells.size(); i++) {
    std::cout << "Cell " << i << std::endl;
    _cells[i].print_particles();
    std::cout << std::endl;
  }
}

Particle Grid::get_first_particle() const { return _first_particle; }

long Grid::get_collisions() {
  long total = 0;
  for (const auto &c : _cells) {
    total += c.collisions;
  }
  return total;
}
