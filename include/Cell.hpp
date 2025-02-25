#ifndef __CELL_HPP__
#define __CELL_HPP__

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

#include "Particle.hpp"
#include <vector>

class Cell {
private:
  std::vector<Particle> _particles;
  std::vector<Particle> _temp_particles;
  double _x, _y;
  long _side;
  double _center_of_mass_x, _center_of_mass_y, _mass;

  void _update_mass();
  void _update_center_of_mass();

public:
  long _collisions = 0;

  Cell(double x, double y, long side);

  void add_particle(Particle &p);
  bool is_particle_inside(Particle &p);
  std::vector<Particle> update_particles(std::vector<Cell> &adjacent_cells, double grid_side);
  void finish_update();
  void print_particles();

  double get_center_of_mass_x();
  double get_center_of_mass_y();
  double get_cell_mass();
  std::vector<Particle>& get_particles();
};

#endif
