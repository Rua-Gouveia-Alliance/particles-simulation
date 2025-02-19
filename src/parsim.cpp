#define _USE_MATH_DEFINES

#include "Particle.hpp"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

unsigned int seed;

void init_r4uni(int input_seed) { seed = input_seed + 987654321; }

double rnd_uniform01() {
  int seed_in = seed;
  seed ^= (seed << 13);
  seed ^= (seed >> 17);
  seed ^= (seed << 5);
  return 0.5 + 0.2328306e-09 * (seed_in + static_cast<int>(seed));
}

double rnd_normal01() {
  double u1, u2, z, result;
  do {
    u1 = rnd_uniform01();
    u2 = rnd_uniform01();
    z = std::sqrt(-2 * std::log(u1)) * std::cos(2 * M_PI * u2);
    result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
  } while (result < 0 || result >= 1);
  return result;
}

void init_particles(long seed, double side, long ncside, long long n_part,
                    std::vector<Particle> &par) {
  double (*rnd01)() = rnd_uniform01;
  long long i;

  if (seed < 0) {
    rnd01 = rnd_normal01;
    seed = -seed;
  }

  init_r4uni(seed);
  par.resize(n_part);

  for (i = 0; i < n_part; i++) {
    par[i]._x = rnd01() * side;
    par[i]._y = rnd01() * side;
    par[i]._vx = (rnd01() - 0.5) * side / ncside / 5.0;
    par[i]._vy = (rnd01() - 0.5) * side / ncside / 5.0;

    par[i]._m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
  }
}

void simulation(std::vector<Particle> &par) {
  // TODO
}

void print_result() {
  // TODO
}

int main(int argc, char *argv[]) {
  double exec_time;
  long long ll;
  std::vector<Particle> particles;

  if (argc != 6) {
    std::cerr << "Usage: " << argv[0]
              << " <seed> <side> <ncside> <n_part> <time_steps>\n";
    return 1;
  }

  try {
    long seed = std::stol(argv[1]);
    double side = std::stod(argv[2]);
    long ncside = std::stol(argv[3]);
    long long n_part = std::stoll(argv[4]);
    long long time_steps = std::stoll(argv[5]);

    init_particles(seed, side, ncside, n_part, particles);

    exec_time = -omp_get_wtime();
    for (ll = 0; ll < time_steps; ll++) {
      simulation(particles);
    }
    exec_time += omp_get_wtime();

    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(); // to stdout
  } catch (const std::exception &e) {
    std::cerr << "Error: Invalid input." << e.what() << "\n";
    return 1;
  }

  return 0;
}
