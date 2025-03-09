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
  #pragma omp parallel for reduction(+:total)
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
  #pragma omp parallel for reduction(+:x_res, y_res)
  for (const auto &p : _particles) {
    x_res += p.m * p.x;
    y_res += p.m * p.y;
  }

  center_of_mass_x = x_res / mass;
  center_of_mass_y = y_res / mass;
}

void Cell::_check_collisions() {
  //TODO otimizar 
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

std::vector<Particle> Cell::update_particles(const std::vector<Cell> &adjacent_cells) {
  std::vector<Particle> new_particles(_particles.size());

  #pragma omp parallel
  {
      std::vector<std::pair<double, double>> local_forces(_particles.size(), {0.0, 0.0});

      #pragma omp for schedule(dynamic) nowait
      for (long long i = 0; i < _particles.size(); i++) {
          const Particle &pi = _particles[i];
          double Gm_i = G * pi.m;

          // Calculate forces within the same cell
          #pragma omp simd
          for (long long j = i + 1; j < _particles.size(); j++) {
              double dx = _particles[j].x - pi.x;
              double dy = _particles[j].y - pi.y;
              double distance_sq = dx * dx + dy * dy;

              double force = Gm_i * _particles[j].m / distance_sq;
              double inv_distance_sqrt = 1.0 / std::sqrt(distance_sq);
              double force_x = force * dx * inv_distance_sqrt;
              double force_y = force * dy * inv_distance_sqrt;

              local_forces[i].first += force_x;
              local_forces[j].first -= force_x;
              local_forces[i].second += force_y;
              local_forces[j].second -= force_y;
          }

          // Calculate forces from adjacent cells
          for (const auto &ac : adjacent_cells) {
              if (ac.mass == 0) continue;

              double dx = ac.center_of_mass_x - pi.x;
              double dy = ac.center_of_mass_y - pi.y;
              double distance_sq = dx * dx + dy * dy;

              double force = Gm_i * ac.mass / distance_sq;
              double inv_distance_sqrt = 1.0 / std::sqrt(distance_sq);

              local_forces[i].first += force * dx * inv_distance_sqrt;
              local_forces[i].second += force * dy * inv_distance_sqrt;
          }
      }

      #pragma omp for schedule(dynamic) nowait
      for (long long i = 0; i < _particles.size(); i++) {
          double ax = local_forces[i].first / _particles[i].m;
          double ay = local_forces[i].second / _particles[i].m;

          new_particles[i] = Particle(_particles[i].m);
          new_particles[i].vx = _particles[i].vx + ax * DELTAT;
          new_particles[i].vy = _particles[i].vy + ay * DELTAT;
          new_particles[i].x = _particles[i].x + _particles[i].vx * DELTAT + 0.5 * ax * DELTAT * DELTAT;
          new_particles[i].y = _particles[i].y + _particles[i].vy * DELTAT + 0.5 * ay * DELTAT * DELTAT;

          if (_particles[i].first_particle)
              new_particles[i].first_particle = true;
      }
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
