#include "Grid.hpp"
#include <iostream>

Grid::Grid(double side, long ncside) : _side(side), _ncside(ncside) {}

void Grid::_add_particle_to_cell(Particle &p) {
  double cell_size = this->_side / this->_ncside;
  long cell_x = static_cast<long>(p._x / cell_size) % this->_ncside;
  long cell_y = static_cast<long>(p._y / cell_size) % this->_ncside;

  if (cell_x < 0)
    cell_x += this->_ncside;
  if (cell_y < 0)
    cell_y += this->_ncside;

  long cell_idx = cell_x + cell_y * this->_ncside;
  this->_cells[cell_idx].add_particle(p);
}

void Grid::add_cell(Cell &c) { this->_cells.push_back(c); }

std::vector<long> Grid::get_adjacent_cells(long ci) {
  std::vector<long> adjacent;
  long cx = ci % this->_ncside;
  long cy = ci / this->_ncside;

  for (long ay = -1; ay <= 1; ay++) {
    for (long ax = -1; ax <= 1; ax++) {
      if (ay == 0 && ax == 0)
        continue;

      long nx = (cx + ax + this->_ncside) % this->_ncside;
      long ny = (cy + ay + this->_ncside) % this->_ncside;

      long idx = nx + ny * this->_ncside;
      adjacent.push_back(idx);
    }
  }

  return adjacent;
}

void Grid::update_cells() {
  for (long i = 0; i < this->_cells.size(); i++) {
    Cell &curr_cell = this->_cells[i];
    double cell_side;

    const std::vector<long> &adjacent_idx = this->get_adjacent_cells(i);
    std::vector<Cell> adjacent_cells;
    for (long j : adjacent_idx) {
      Cell new_cell = this->_cells[j];
      cell_side = curr_cell._side;

      // wrapping in x direction
      if (new_cell._x >= curr_cell._x + cell_side * 2) {
        new_cell._center_of_mass_x -= this->_side;
      } else if (new_cell._x < curr_cell._x - cell_side) {
        new_cell._center_of_mass_x += this->_side;
      }

      // wrapping in y direction
      if (new_cell._y >= curr_cell._y + cell_side * 2) {
        new_cell._center_of_mass_y -= this->_side;
      } else if (new_cell._y < curr_cell._y - cell_side) {
        new_cell._center_of_mass_y += this->_side;
      }

      adjacent_cells.push_back(new_cell);
    }

    std::vector<Particle> new_particles =
        curr_cell.update_particles(adjacent_cells);

    for (auto &p : new_particles) {
      // small optimization
      if (curr_cell.is_particle_inside(p))
        curr_cell.add_particle(p);
      else
        this->_add_particle_to_cell(p);
    }
  }

  for (auto &c : this->_cells) {
    c.finish_update();
  }
}

void Grid::print_cells() {
  for (long i = 0; i < this->_cells.size(); i++) {
    std::cout << "Cell " << i << std::endl;
    this->_cells[i].print_particles();
    std::cout << std::endl;
  }
}

Particle Grid::get_first_particle() {
  // TODO this is not efficient at all, maybe improve
  for (auto &c : this->_cells) {
    for (auto &p : c.get_particles()) {
      if (p._first_particle)
        return p;
    }
  }
  return Particle();
}

long Grid::get_collisions() {
  long total = 0;
  for (auto &c : this->_cells) {
    total += c._collisions;
  }
  return total;
}
