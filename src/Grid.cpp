#include "Grid.hpp"

Grid::Grid(long side, long ncside) : _side(side), _ncside(ncside) {}

void Grid::add_cell(Cell c) { this->_cells.push_back(c); }

std::vector<long> Grid::get_adjacent_cells(long ci) {
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

    this->_cells[i].update_particles(adjacent_cells);
  }

  for (auto &c : this->_cells) {
    c.finish_update();
  }
}
