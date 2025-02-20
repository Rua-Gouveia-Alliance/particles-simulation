#ifndef __GRID_HPP__
#define __GRID_HPP__

#include "Cell.hpp"
#include <vector>

class Grid {
private:
  std::vector<Cell> _cells;

public:
  Grid();

  void add_cell(Cell c);
  void update_cells();
};

#endif
