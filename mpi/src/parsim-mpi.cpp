#include "Cell.hpp"
#include "Grid.hpp"
#include "Particle.hpp"
#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mpi.h>
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
void compute_grid_partition(long ncside, int size, int rank, long &start_row,
                            long &end_row) {
  long rows_per_proc = ncside / size;
  long extra_rows = ncside % size;
  start_row =
      rank * rows_per_proc + std::min(static_cast<long>(rank), extra_rows);
  end_row = start_row + rows_per_proc + (rank < extra_rows ? 1 : 0);
}

void init_particles(long seed, double side, long ncside, long long n_part,
                    std::vector<Particle> &local_par, int rank, int size) {
  double (*rnd01)() = rnd_uniform01;
  long long i;

  if (seed < 0) {
    rnd01 = rnd_normal01;
    seed = -seed;
  }

  init_r4uni(seed);
  local_par.clear();

  long start_row, end_row;
  compute_grid_partition(ncside, size, rank, start_row, end_row);

  /*for (i = 0; i < n_part; i++) {
    par[i]._x = rnd01() * side;
    par[i]._y = rnd01() * side;
    par[i]._vx = (rnd01() - 0.5) * side / ncside / 5.0;
    par[i]._vy = (rnd01() - 0.5) * side / ncside / 5.0;

    par[i]._m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
  }
  par[0]._first_particle = true;*/

  double cell_size = side / ncside;

  double y_min = start_row * cell_size;
  double y_max = end_row * cell_size;
  for (i = 0; i < n_part; i++) {
    double x = rnd01() * side;
    double y = rnd01() * side;
    if (y >= y_min && y < y_max) {
      Particle p;
      p._x = x;
      p._y = y;
      p._vx = (rnd01() - 0.5) * side / ncside / 5.0;
      p._vy = (rnd01() - 0.5) * side / ncside / 5.0;
      p._m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
      if (i == 0) {
        p._first_particle = true;
      }
      local_par.push_back(p);
    }
  }
}

Grid init_grid(double side, long ncside, std::vector<Particle> &pv, int rank,
               int size) {

  long start_row, end_row;
  compute_grid_partition(ncside, size, rank, start_row, end_row);

  double cell_size = side / ncside;

  Grid local_grid(side, ncside);

  // create cells
  /* for (long i = 0; i < ncside; i++) {
     double y = i * cell_size;
     for (long j = 0; j < ncside; j++) {
       double x = j * cell_size;

       Cell cell(x, y, cell_size, rank);
       grid.add_cell(cell);
     }
   }*/
  for (long i = start_row; i < end_row; i++) {
    double y = i * cell_size;
    for (long j = 0; j < ncside; j++) {
      double x = j * cell_size;

      Cell cell(x, y, cell_size, rank);
      local_grid.add_cell(cell);
    }
  }

  // assign particles to corresponding cells
  for (auto &p : pv) {
    long cell_x = static_cast<long>(p._x / cell_size);
    long cell_y = static_cast<long>(p._y / cell_size);

    if (cell_y >= start_row && cell_y < end_row) {
      long cell_idx = cell_x + (cell_y - start_row) * ncside;
      local_grid._cells[cell_idx].add_particle(p);
    }
  }

  // initialize cell
  for (auto &c : local_grid._cells) {
    c.finish_update();
  }

  return local_grid;
}

void simulation(Grid &grid, long long time_steps) {
  for (long long ll = 0; ll < time_steps; ll++) {
    // std::cout << "Round " << ll << std::endl; // debug
    grid.update_cells();
    // grid.print_cells(); // debug
  }
}

void print_result(Grid &g) {
  // g.print_cells(); // debug
  Particle pf = g.get_first_particle();
  fprintf(stdout, "%.3f %.3f\n", pf._x, pf._y);

  std::cout << g.get_collisions() << std::endl;
}

int main(int argc, char *argv[]) {
  double exec_time;
  MPI_Init(&argc, &argv); // init

  std::vector<Particle> local_par;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    init_particles(seed, side, ncside, n_part, local_par, rank, size);

    // debug
    // std::cout << "INITIAL GRID" << std::endl;
    // grid.print_cells();

    exec_time = -omp_get_wtime();
    Grid grid = init_grid(side, ncside, local_par, rank, size);
    simulation(grid, time_steps);
    exec_time += omp_get_wtime();

    fprintf(stderr, "%.1fs\n", exec_time);
    print_result(grid); // to stdout
  } catch (const std::exception &e) {
    std::cerr << "Error: Invalid input." << e.what() << "\n";
    MPI_Finalize();
    return 1;
  }
  MPI_Finalize();
  return 0;
}
