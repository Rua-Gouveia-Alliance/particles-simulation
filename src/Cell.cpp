#include "Cell.hpp"
#include "Particle.hpp"

#include <cmath>
#include <iostream>
#include <vector>

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

  if (this->_mass == 0) {
    this->_center_of_mass_x = -1;
    this->_center_of_mass_y = -1;
    return;
  }

  for (Particle &p : this->_particles) {
    x_res += p._m * p._x;
    y_res += p._m * p._y;
  }

  this->_center_of_mass_x = x_res / this->_mass;
  this->_center_of_mass_y = y_res / this->_mass;
}

void Cell::add_particle(Particle &p) { this->_temp_particles.push_back(p); }

bool Cell::is_particle_inside(Particle &p) {
  return this->_x <= p._x && p._x <= this->_x + this->_side &&
         this->_y <= p._y && p._y <= this->_y + this->_side;
}

std::vector<Particle>
Cell::update_particles(std::vector<Cell> &adjacent_cells) {
  std::vector<Particle> new_particles;

  // 2 vectors of force to save and only calculate A->B and avoid B->A
  std::vector<std::pair<double, double>> forces(this->_particles.size(),
                                                {0.0, 0.0});

  for (long long i = 0; i < this->_particles.size(); i++) {
    Particle &pi = this->_particles[i], new_particle;
    double force;
    double dx, dy, distance_sq, distance_sqrt, inv_distance_sqrt;
    double force_x = 0.0, force_y = 0.0;
    double ax, ay;
    bool collided = false;
    double Gm_i = G * this->_particles[i]._m;

    // calculate resulting force for particles inside same cell
    for (long long j = i + 1; j < this->_particles.size(); j++) {

      dx = pi._x - this->_particles[j]._x;
      dy = pi._y - this->_particles[j]._y;
      distance_sq = dx * dx + dy * dy;

      if (distance_sq < EPSILON2) {
        collided = true;
        this->_collisions++;
        break;
      }

      force = Gm_i * this->_particles[j]._m;
      force /= distance_sq;

      inv_distance_sqrt = 1.0 / std::sqrt(distance_sq);
      force_x += force * (dx * inv_distance_sqrt);
      force_y += force * (dy * inv_distance_sqrt);

      forces[i].first += force_x;
      forces[j].first -= force_x;
      forces[i].second += force_y;
      forces[j].second -= force_y;
    }

    if (collided)
      continue;

    // calculate resulting force for adjacent cells
    for (auto &ac : adjacent_cells) {
      if (ac._mass == 0)
        continue;

      dx = pi._x - ac._center_of_mass_x;
      dy = pi._y - ac._center_of_mass_y;

      distance_sq = dx * dx + dy * dy;

      force = Gm_i * ac._mass;
      force /= distance_sq;

      inv_distance_sqrt = 1.0 / std::sqrt(distance_sq);
      forces[i].first += force * (dx * inv_distance_sqrt);
      forces[i].second += force * (dy * inv_distance_sqrt);
    }

    // calculate new acceleration, velocity, position
    new_particle = Particle(pi._m);

    ax = forces[i].first / pi._m;
    ay = forces[i].second / pi._m;

    new_particle._vx = pi._vx + ax * DELTAT;
    new_particle._vy = pi._vy + ay * DELTAT;

    new_particle._x = pi._x + pi._vx * DELTAT + 0.5 * ax * (DELTAT * DELTAT);
    new_particle._y = pi._y + pi._vy * DELTAT + 0.5 * ay * (DELTAT * DELTAT);

    if (pi._first_particle)
      new_particle._first_particle = true;

    new_particles.push_back(new_particle);
  }

  return new_particles;
}

void Cell::finish_update() {
  this->_particles = this->_temp_particles;
  this->_temp_particles = std::vector<Particle>();
  this->_update_mass();
  this->_update_center_of_mass();
}

void Cell::print_particles() {
  for (auto &p : this->_particles) {
    p.print_info();
  }
  std::cout << "Center of Mass: " << this->_center_of_mass_x << ", "
            << this->get_center_of_mass_y() << std::endl;
  std::cout << "Mass: " << this->_mass << std::endl;
}

double Cell::get_center_of_mass_x() { return this->_center_of_mass_x; }
double Cell::get_center_of_mass_y() { return this->_center_of_mass_y; }
double Cell::get_cell_mass() { return this->_mass; }
std::vector<Particle> &Cell::get_particles() { return this->_particles; }
