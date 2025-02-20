#include "Grid.hpp"

Grid::Grid(long side, long ncside) : _side(side), _ncside(ncside) {}

void Grid::_add_particle_to_cell(Particle &p) {
  double cell_size = this->_side / this->_ncside;
  long cell_x = static_cast<long>(p._x / cell_size);
  long cell_y = static_cast<long>(p._y / cell_size);

  long cell_idx = cell_x + cell_y * this->_ncside;
  // TODO wrap around
  this->_cells[cell_idx].add_particle(p);
}

void Grid::add_cell(Cell &c) { this->_cells.push_back(c); }

std::vector<long> Grid::get_adjacent_cells(long ci) {
  // TODO wrap around
  std::vector<long> adjacent;
  long cx = ci % this->_ncside;
  long cy = ci / this->_ncside;

  for (long ay = -1; ay <= 1; ay++) {
    for (long ax = -1; ax <= 1; ax++) {
      if (ay == 0 && ax == 0)
        continue;

      long nx = cx + ax;
      long ny = cy + ay;

      if (ny >= 0 && ny < this->_ncside && nx >= 0 && nx < this->_ncside) {
        long idx = nx + ny * this->_ncside;
        adjacent.push_back(idx);
      }
    }
  }

  return adjacent;
}

void Grid::update_cells() {
  for (long i = 0; i < this->_cells.size(); i++) {
    std::vector<long> adjacent_idx = this->get_adjacent_cells(i);
    std::vector<Cell> adjacent_cells;
    for (long j : adjacent_idx) {
      adjacent_cells.push_back(this->_cells[j]);
    }

    std::vector<Particle> new_particles =
        this->_cells[i].update_particles(adjacent_cells);

    for (auto &p : new_particles) {
      // TODO i dont know if this check improves performance or if its worse
      if (this->_cells[i].is_particle_inside(p))
        this->_cells[i].add_particle(p);
      else
        this->_add_particle_to_cell(p);
    }
  }

  for (auto &c : this->_cells) {
    c.finish_update();
  }
}
