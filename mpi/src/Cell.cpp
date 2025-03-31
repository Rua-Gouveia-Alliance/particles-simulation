#include "Cell.hpp"
#include "Mass.hpp"
#include "Particle.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

void Cell::update_mass() {
  double total = 0;

  for (const auto &p : _particles) {
    total += p.m;
  }
  mass.val = total;

  double x_res = 0;
  double y_res = 0;

  if (mass.val == 0) {
    mass.x = -1;
    mass.y = -1;
    return;
  }

  for (const auto &p : _particles) {
    x_res += p.m * p.x;
    y_res += p.m * p.y;
  }

  mass.x = x_res / mass.val;
  mass.y = y_res / mass.val;
}

void Cell::check_collisions() {
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

      dx = pi.x - pj.x;
      dy = pi.y - pj.y;
      distance_sq = dx * dx + dy * dy;

      if (distance_sq < EPSILON2) {
        if (!collided[i] && !collided[j])
          collisions++;
        collided[i] = true;
        collided[j] = true;
      }
    }

    if (!collided[i])
      final_particles.push_back(_particles[i]);
  }

  _particles = final_particles;
}

std::vector<Particle>
Cell::update_particles(const std::vector<Mass> &adjacent_masses) {
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

      if (pi.id == 0) {
        std::cout << "Force by particle " << _particles[j].id << ", " << force_x
                  << ", " << force_y << "\n";
      }

      forces[i].first += force_x;
      forces[j].first -= force_x;
      forces[i].second += force_y;
      forces[j].second -= force_y;
    }

    // calculate resulting force for adjacent cells
    for (const auto &mass : adjacent_masses) {
      if (mass.val == 0)
        continue;

      dx = mass.x - pi.x;
      dy = mass.y - pi.y;
      distance_sq = dx * dx + dy * dy;

      force = Gm_i * mass.val;
      force /= distance_sq;

      inv_distance_sqrt = 1.0 / std::sqrt(distance_sq);
      forces[i].first += force * (dx * inv_distance_sqrt);
      forces[i].second += force * (dy * inv_distance_sqrt);

      if (pi.id == 0)
        std::cout << "Force by cell " << mass.id << ", "
                  << force * (dx * inv_distance_sqrt) << ", "
                  << force * (dy * inv_distance_sqrt) << "\n";
    }

    // calculate new acceleration, velocity, position
    new_particle = {.m = pi.m};

    ax = forces[i].first / pi.m;
    ay = forces[i].second / pi.m;

    new_particle.vx = pi.vx + ax * DELTAT;
    new_particle.vy = pi.vy + ay * DELTAT;

    new_particle.x = pi.x + pi.vx * DELTAT + 0.5 * ax * (DELTAT * DELTAT);
    new_particle.y = pi.y + pi.vy * DELTAT + 0.5 * ay * (DELTAT * DELTAT);

    if (pi.first_particle)
      new_particle.first_particle = true;

    new_particle.id = pi.id;
    new_particles.push_back(new_particle);
  }

  return new_particles;
}

void Cell::finish_update() {
  _particles = _temp_particles;
  _temp_particles = std::vector<Particle>();
}
