#include "Cell.hpp"
#include "Grid.hpp"
#define _USE_MATH_DEFINES

#include "Grid.hpp"
#include "Particle.hpp"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

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

Grid init_grid(double side, long ncside, std::vector<Particle> &pv) {
  Grid grid(side, ncside);
  double cell_size = side / ncside;

  // create cells
  for (long i = 0; i < ncside; i++) {
    for (long j = 0; j < ncside; j++) {
      double x = i * cell_size;
      double y = j * cell_size;

      Cell cell(x, y, cell_size);
      grid.add_cell(cell);
    }
  }

  // assign particles to corresponding cells
  for (auto &p : pv) {
    long cell_x = static_cast<long>(p._x / cell_size);
    long cell_y = static_cast<long>(p._y / cell_size);

    long cell_idx = cell_x + cell_y * ncside;
    grid._cells[cell_idx].add_particle(p);
  }

  // initialize cell
  for (auto &c : grid._cells) {
    c.finish_update();
  }

  return grid;
}

void simulation(Grid &grid, long long time_steps) {
  for (long long ll = 0; ll < time_steps; ll++) {
    grid.update_cells();
  }
}

void print_result(Grid &g) {
  // TODO this should print coords of particle 0 (two real values using three
  // decimal digits) in the first line, and the number of particles that
  // collided in the second line
  g.print_cells();
}

int main(int argc, char *argv[]) {
  double exec_time;
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

    Grid grid = init_grid(side, ncside, particles);
    std::cout << "INITIAL GRID" << std::endl;
    grid.print_cells();

    exec_time = -omp_get_wtime();
    simulation(grid, time_steps);
    exec_time += omp_get_wtime();

    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(grid); // to stdout
  } catch (const std::exception &e) {
    std::cerr << "Error: Invalid input." << e.what() << "\n";
    return 1;
  }

  return 0;
}
