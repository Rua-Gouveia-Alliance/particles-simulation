#include "Grid.hpp"
#include <iostream>

Grid::Grid(double side, long ncside)
    : _side(side), _ncside(ncside), _first_particle(_default_first_particle) {}

void Grid::_add_particle_to_cell(Particle &p) {
  double temp_px, temp_py;

  // wrap particle around
  do {
    temp_px = p._x;
    temp_py = p._y;

    if (p._x < 0)
      p._x += _side;
    else if (p._x >= _side)
      p._x -= _side;

    if (p._y < 0)
      p._y += _side;
    else if (p._y >= _side)
      p._y -= _side;
  } while (temp_px != p._x || temp_py != p._y);

  // calculate cell index
  double cell_size = _side / _ncside;
  long cell_x = static_cast<long>(p._x / cell_size) % _ncside;
  long cell_y = static_cast<long>(p._y / cell_size) % _ncside;

  long cell_idx = cell_x + cell_y * _ncside;
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
  for (long i = 0; i < _cells.size(); i++) {
    Cell &curr_cell = _cells[i];
    double cell_side;

    const std::vector<long> &adjacent_idx = get_adjacent_cells(i);
    std::vector<Cell> adjacent_cells;
    for (long j : adjacent_idx) {
      Cell new_cell = _cells[j];
      cell_side = curr_cell._side;

      // wrapping in x direction
      if (new_cell._x >= curr_cell._x + cell_side * 2) {
        new_cell._center_of_mass_x -= _side;
      } else if (new_cell._x < curr_cell._x - cell_side) {
        new_cell._center_of_mass_x += _side;
      }

      // wrapping in y direction
      if (new_cell._y >= curr_cell._y + cell_side * 2) {
        new_cell._center_of_mass_y -= _side;
      } else if (new_cell._y < curr_cell._y - cell_side) {
        new_cell._center_of_mass_y += _side;
      }

      adjacent_cells.push_back(new_cell);
    }

    std::vector<Particle> new_particles =
        curr_cell.update_particles(adjacent_cells);

    for (auto &p : new_particles) {
      if (p._first_particle)
        _first_particle = p;

      _add_particle_to_cell(p);
    }
  }

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
    total += c._collisions;
  }
  return total;
}
