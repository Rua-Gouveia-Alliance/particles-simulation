#include "Cell.hpp"
#include "Particle.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <vector>

Cell::Cell(double x, double y, double side) : x(x), y(y), side(side){};

void Cell::_update_mass() {
  double total = 0;
  for (const auto &p : _particles) {
    total += p.m;
  }
  mass = total;
}

void Cell::_update_center_of_mass() {
  double x_res = 0;
  double y_res = 0;

  if (mass == 0) {
    center_of_mass_x = -1;
    center_of_mass_y = -1;
    return;
  }

  for (const auto &p : _particles) {
    x_res += p.m * p.x;
    y_res += p.m * p.y;
  }

  center_of_mass_x = x_res / mass;
  center_of_mass_y = y_res / mass;
}

void Cell::_check_collisions() {
  std::vector<Particle> final_particles;
  long long p_count = _particles.size();
  std::vector<bool> collided(p_count, false);
  double dx, dy, distance_sq, distance_sq_x_k, distance_sq_j_k;
  long long i;

  for (i = 0; i < p_count; i++) {
    if(collided[i]) continue;
    bool found_colision = false;

    const Particle &pi = _particles[i];

    for (long long j = i + 1; j < p_count; j++) {
      if(collided[j]) continue;

      const Particle &pj = _particles[j];

      dx = pi.x - pj.x;
      dy = pi.y - pj.y;
      distance_sq = dx * dx + dy * dy;

      if (distance_sq < EPSILON2) {
        bool three_particles_collided = false;
        for(long long k = j + 1; k < p_count; k++)
        {
          if (collided[k]) continue;
          const Particle &pk = _particles[k];

          dx = pi.x - pk.x;
          dy = pi.y - pk.y;

          distance_sq_x_k = dx * dx + dy * dy;

          dx = pj.x - pk.x;
          dy = pj.y - pk.y;

          distance_sq_j_k =  dx * dx + dy * dy;

          if(distance_sq_x_k < EPSILON2 && distance_sq_j_k < EPSILON2)
          {
            // 3 particle collided
            collisions++;
            collided[i] = true;
            collided[j] = true;
            collided[k] = true;
            found_colision = true;
            break;
          }
         
      }
      if(!found_colision)
      {
        //only 2 collided 
        collided[i] = true;
        collided[j] = true;
        collisions ++;
        found_colision = true;
      }
      break;
      }
    }
    
    if (!collided[i])
      final_particles.push_back(_particles[i]);
  }

  _particles = final_particles;
}

void Cell::add_particle(Particle &p) { _temp_particles.push_back(p); }

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
    Gm_i = G * _particles[i].m;

    // calculate resulting force for particles inside same cell
    for (long long j = i + 1; j < _particles.size(); j++) {
      dx = _particles[j].x - pi.x;
      dy = _particles[j].y - pi.y;
      distance_sq = dx * dx + dy * dy;

      force = Gm_i * _particles[j].m;
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
    for (const auto &ac : adjacent_cells) {
      if (ac.mass == 0)
        continue;

      dx = ac.center_of_mass_x - pi.x;
      dy = ac.center_of_mass_y - pi.y;
      distance_sq = dx * dx + dy * dy;

      force = Gm_i * ac.mass;
      force /= distance_sq;

      inv_distance_sqrt = 1.0 / std::sqrt(distance_sq);
      forces[i].first += force * (dx * inv_distance_sqrt);
      forces[i].second += force * (dy * inv_distance_sqrt);
    }

    // calculate new acceleration, velocity, position
    new_particle = Particle(pi.m);

    ax = forces[i].first / pi.m;
    ay = forces[i].second / pi.m;

    new_particle.vx = pi.vx + ax * DELTAT;
    new_particle.vy = pi.vy + ay * DELTAT;

    new_particle.x = pi.x + pi.vx * DELTAT + 0.5 * ax * (DELTAT * DELTAT);
    new_particle.y = pi.y + pi.vy * DELTAT + 0.5 * ay * (DELTAT * DELTAT);

    if (pi.first_particle)
      new_particle.first_particle = true;

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

void Cell::print_particles() const {
  for (const auto &p : _particles) {
    p.print_info();
  }
  std::cout << "Center of Mass: " << center_of_mass_x << ", "
            << center_of_mass_y << std::endl;
  std::cout << "Mass: " << mass << std::endl;
}

const std::vector<Particle> &Cell::get_particles() { return _particles; }
