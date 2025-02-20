#include "Cell.hpp"

Cell::Cell(double x, double y, long side) : _x(x), _y(y), _side(side){};

void Cell::_update_center_of_mass() {
  double x_res = 0;
  double y_res = 0;

  for (Particle &p : this->_particles) {
    x_res += p._m * p._x;
    y_res += p._m * p._y;
  }

  this->_center_of_mass_x = x_res / this->_mass;
  this->_center_of_mass_y = y_res / this->_mass;
}

void Cell::_remove_collided() {
  // TODO
}

void Cell::add_particle(Particle &p) { this->_particles.push_back(p); }

bool Cell::is_particle_inside(Particle &p) {
  double bb_ix = this->_x, bb_iy = this->_y, bb_ax = this->_x + this->_side,
         bb_ay = this->_y + this->_side, p_x = p._x, p_y = p._y;

  if (bb_ix <= p_x && p_x <= bb_ax && bb_iy <= p_y && p_y <= bb_ay) {
    return true;
  }
  return false;
}

void update_particles(std::vector<Cell> &adjacent_cells) {
  // TODO
}

void Cell::update_collisions() {
  // TODO
}

void Cell::finish_update() {
  this->_update_center_of_mass();
  this->_remove_collided();
}

double Cell::get_center_of_mass_x() { return this->_center_of_mass_x; }
double Cell::get_center_of_mass_y() { return this->_center_of_mass_y; }
double Cell::get_cell_mass() { return this->_mass; }
