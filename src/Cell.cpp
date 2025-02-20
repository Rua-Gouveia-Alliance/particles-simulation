#include "Cell.hpp"

#include <cmath>

Cell::Cell(double x, double y, long side) : _x(x), _y(y), _side(side){};

void Cell::_update_mass() {
  double total = 0;
  for (Particle &p : this->_particles) {
    total += p._m;
  }
  this->_mass = total;
}

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
  for (long long i = 0; i < this->_particles.size(); i++) {
    if (this->_particles[i]._collided) {
      this->_particles.erase(this->_particles.begin() + i);
      i--;
    }
  }
}

void Cell::add_particle(Particle &p) { this->_temp_particles.push_back(p); }

bool Cell::is_particle_inside(Particle &p) {
  double bb_ix = this->_x, bb_iy = this->_y, bb_ax = this->_x + this->_side,
         bb_ay = this->_y + this->_side, p_x = p._x, p_y = p._y;

  if (bb_ix <= p_x && p_x <= bb_ax && bb_iy <= p_y && p_y <= bb_ay) {
    return true;
  }
  return false;
}

void Cell::update_particles(std::vector<Cell> &adjacent_cells) {
  for (long long i = 0; i < this->_particles.size(); i++) {
    Particle pi = this->_particles[i];
    long long force = 0;

    // TODO calculate Fx and Fy, not only total F
    // TODO check for collisions between particles
    // calculate resulting force for particles inside same cell
    for (long long j = 0; j < this->_particles.size(); j++) {
      if (i == j)
        continue;

      // TODO can be improved, only calculate the force for A,B and B,A once
      force += G * this->_particles[i]._m * this->_particles[j]._m;
      force /= std::pow(pi._x - this->_particles[j]._x, 2) +
               std::pow(pi._y - this->_particles[j]._y, 2);
    }

    // calculate resulting force for adjacent cells
    for (auto &ac : adjacent_cells) {
      force += G * this->_particles[i]._m * ac._mass;
      force /= std::pow(pi._x - ac._center_of_mass_x, 2) +
               std::pow(pi._y - ac._center_of_mass_y, 2);
    }

    // TODO calculate new acceleration, velocity, position

    // TODO create new particle with new acc, vel, pos and add it to this Cell's
    // temp particles, or another cell if it has changed
  }
}

void Cell::finish_update() {
  this->_particles = this->_temp_particles;
  this->_remove_collided();
  this->_update_mass();
  this->_update_center_of_mass();
}

double Cell::get_center_of_mass_x() { return this->_center_of_mass_x; }
double Cell::get_center_of_mass_y() { return this->_center_of_mass_y; }
double Cell::get_cell_mass() { return this->_mass; }
