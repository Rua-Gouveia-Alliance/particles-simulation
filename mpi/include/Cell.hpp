#ifndef __CELL_HPP__
#define __CELL_HPP__

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

#include "Mass.hpp"
#include "Particle.hpp"
#include <vector>

class Cell {
private:
  std::vector<Particle> _particles;
  std::vector<Particle> _temp_particles;

  void _update_mass();
  void _update_center_of_mass();
  void _check_collisions();

public:
  Mass _mass;
  long _collisions = 0;
  double _x, _y, _side;

  Cell(long id, double x, double y, double side);

  void add_particle(Particle &p);
  std::vector<Particle>
  update_particles(const std::vector<Cell> &adjacent_cells);
  void finish_update();
  void print_particles() const;

  long id() const;
  double get_center_of_mass_x() const;
  double get_center_of_mass_y() const;
  double get_cell_mass() const;
  const std::vector<Particle> &get_particles();
};

#endif
