#include "Cell.hpp"
#include "Particle.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <vector>

Cell::Cell(double x, double y, double side) : _x(x), _y(y), _side(side){};

void Cell::_update_mass() {
  double total = 0;
  for (const auto &p : _particles) {
    total += p._m;
  }
  _mass = total;
}

void Cell::_update_center_of_mass() {
  double x_res = 0;
  double y_res = 0;

  if (_mass == 0) {
    _center_of_mass_x = -1;
    _center_of_mass_y = -1;
    return;
  }

  for (const auto &p : _particles) {
    x_res += p._m * p._x;
    y_res += p._m * p._y;
  }

  _center_of_mass_x = x_res / _mass;
  _center_of_mass_y = y_res / _mass;
}

void Cell::_check_collisions() {
  std::vector<Particle> final_particles;
  long long p_count = _particles.size();
  std::vector<bool> collided(p_count, false);
  double dx, dy, distance_sq;
  long long i;

  for (i = 0; i < p_count; i++) {
    const Particle &pi = _particles[i];

    for (long long j = i + 1; j < p_count; j++) {
      if (i == j)
        continue;

      const Particle &pj = _particles[j];

      dx = pi._x - pj._x;
      dy = pi._y - pj._y;
      distance_sq = dx * dx + dy * dy;

      if (distance_sq < EPSILON2) {
        // debug
        // std::cout << "COLLIDED" << std::endl;
        // pi.print_info();
        // pj.print_info();
        // std::cout << std::endl;

        _collisions++;
        collided[i] = true;
        collided[j] = true;
      }
    }

    if (!collided[i])
      final_particles.push_back(_particles[i]);
  }

  _particles = final_particles;
}

void Cell::add_particle(Particle &p) { _temp_particles.push_back(p); }

bool Cell::is_particle_inside(Particle &p) {
  return _x <= p._x && p._x <= _x + _side && _y <= p._y && p._y <= _y + _side;
}

std::vector<Particle>
Cell::update_particles(const std::vector<Cell> &adjacent_cells) {
  Particle new_particle;
  std::vector<Particle> new_particles;
  std::vector<std::pair<double, double>> forces(_particles.size(), {0.0, 0.0});
  double force;
  double dx, dy, distance_sq, inv_distance_sqrt;
  double force_x, force_y;
  double ax, ay;
  double Gm_i;

  for (long long i = 0; i < _particles.size(); i++) {
    const Particle &pi = _particles[i];
    Gm_i = G * _particles[i]._m;

    // calculate resulting force for particles inside same cell
    for (long long j = i + 1; j < _particles.size(); j++) {
      dx = _particles[j]._x - pi._x;
      dy = _particles[j]._y - pi._y;
      distance_sq = dx * dx + dy * dy;

      force = Gm_i * _particles[j]._m;
      force /= distance_sq;

      inv_distance_sqrt = 1.0 / std::sqrt(distance_sq);
      force_x = force * (dx * inv_distance_sqrt);
      force_y = force * (dy * inv_distance_sqrt);

      forces[i].first += force_x;
      forces[j].first -= force_x;
      forces[i].second += force_y;
      forces[j].second -= force_y;
    }

    // calculate resulting force for adjacent cells
    for (auto &ac : adjacent_cells) {
      if (ac._mass == 0)
        continue;

      dx = ac._center_of_mass_x - pi._x;
      dy = ac._center_of_mass_y - pi._y;
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
  _particles = _temp_particles;
  _temp_particles = std::vector<Particle>();
  _check_collisions();
  _update_mass();
  _update_center_of_mass();
}

void Cell::print_particles() {
  for (const auto &p : _particles) {
    p.print_info();
  }
  std::cout << "Center of Mass: " << _center_of_mass_x << ", "
            << get_center_of_mass_y() << std::endl;
  std::cout << "Mass: " << _mass << std::endl;
}

double Cell::get_center_of_mass_x() { return _center_of_mass_x; }
double Cell::get_center_of_mass_y() { return _center_of_mass_y; }
double Cell::get_cell_mass() { return _mass; }
std::vector<Particle> &Cell::get_particles() { return _particles; }
