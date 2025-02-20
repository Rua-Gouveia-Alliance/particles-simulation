#ifndef __CELL_HPP__
#define __CELL_HPP__

#include "Particle.hpp"
#include <vector>

class Cell {
private:
  std::vector<Particle> _particles;
  double _x, _y;
  long side;
  double center_of_mass;

  void update_center_of_mass();
  void remove_collided();

public:
  Cell(double x, double y, long side);

  void add_particle(Particle p);
  bool is_particle_inside(Particle p);
  void update_collisions();
  void finish_update();

  double get_center_of_mass();
  double get_cell_mass();
};

#endif
