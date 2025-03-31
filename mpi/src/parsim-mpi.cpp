#include "Cell.hpp"
#include "CellLocation.hpp"
#include "PartialCell.hpp"
#include "PartialGrid.hpp"
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

void init_particles(long seed, double side, long ncside, long long n_part,
                    std::vector<Particle> &par, int rank) {
  double (*rnd01)() = rnd_uniform01;
  long long i;

  if (seed < 0) {
    rnd01 = rnd_normal01;
    seed = -seed;
  }

  init_r4uni(seed);
  par.resize(n_part);

  for (i = 0; i < n_part; i++) {
    par[i].x = rnd01() * side;
    par[i].y = rnd01() * side;
    par[i].vx = (rnd01() - 0.5) * side / ncside / 5.0;
    par[i].vy = (rnd01() - 0.5) * side / ncside / 5.0;

    par[i].m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
    par[i].id = i;

    if (rank == 0) {
      fprintf(stdout, "Particle %lli: ", i);
      fprintf(stdout, "mass=%#.6f ", par[i].m);
      fprintf(stdout, "x=%#.6f ", par[i].x);
      fprintf(stdout, "y=%#.6f ", par[i].y);
      fprintf(stdout, "vx=%#.6f ", par[i].vx);
      fprintf(stdout, "vy=%#.6f\n", par[i].vy);
    }
  }
  par[0].first_particle = true;
}

std::vector<CellLocation> partition_grid(int nprocs, long ncside,
                                         double cell_size) {
  std::vector<CellLocation> partition;
  partition.reserve(ncside * ncside);

  for (long i = 0; i < ncside; ++i) {
    double y = i * cell_size;
    for (long j = 0; j < ncside; ++j) {
      double x = j * cell_size;
      partition.push_back(CellLocation(x, y, i % nprocs));
    }
  }

  return partition;
}

// TODO: this can probably be optimized
PartialGrid init_grid(int rank, int nprocs, double side, long ncside,
                      std::vector<Particle> &pv) {
  PartialGrid grid(rank, ncside - 1, side, ncside);
  double cell_size = side / ncside;
  std::vector<CellLocation> partition =
      partition_grid(nprocs, ncside, cell_size);

  std::vector<bool> grid_ranks = std::vector<bool>(nprocs, false);
  for (long i = 0; i < partition.size(); ++i) {
    std::vector<int> adj;
    bool owner_is_adjacent = false;
    CellLocation &loc = partition[i];
    std::vector<bool> ranks = std::vector<bool>(nprocs, false);
    std::vector<long> adj_cells = PartialGrid::get_adjacent_cells(i, ncside);

    adj.reserve(adj_cells.size());
    for (const auto &id : adj_cells) {
      if (partition[id].rank == rank) {
        owner_is_adjacent = true;
        continue;
      }

      if (!ranks[partition[id].rank]) {
        ranks[partition[id].rank] = true;
        adj.push_back(partition[id].rank);
      }
      if (!grid_ranks[partition[id].rank]) {
        grid_ranks[partition[id].rank] = true;
        grid.add_adjacent_rank(partition[id].rank);
      }
      if (!partition[id].counted) {
        partition[id].counted = true;
        grid.increment_adjacent_rank(partition[id].rank);
      }
    }

    if (loc.rank == rank) {
      grid.add_local_cell(Cell(i, adj, loc.x, loc.y, cell_size));
    } else if (owner_is_adjacent) {
      grid.add_adjacent_cell(PartialCell(i, loc.rank, loc.x, loc.y, cell_size));
    }
  }

  // assign particles to corresponding cells
  for (auto &p : pv) {
    long cell_x = static_cast<long>(p.x / cell_size);
    long cell_y = static_cast<long>(p.y / cell_size);
    long cell_idx = cell_x + cell_y * ncside;
    if (partition[cell_idx].rank == rank)
      grid.local_cells.at(cell_idx).add_particle(p);
  }

  // initialize cell
  for (auto &it : grid.local_cells) {
    it.second.finish_update();
    it.second.check_collisions();
  }

  return grid;
}

void print_result(PartialGrid &g) {
  Particle pf = g.get_first_particle();
  fprintf(stdout, "%.3f %.3f\n", pf.x, pf.y);

  std::cout << g.get_collisions() << std::endl;
}

int main(int argc, char *argv[]) {
  double exec_time;
  std::vector<Particle> particles;
  int rank, nprocs;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc != 6) {
    if (rank == 0) {
      std::cerr << "Usage: " << argv[0]
                << " <seed> <side> <ncside> <n_part> <time_steps>\n";
    }
    MPI_Finalize();
    return 1;
  }

  try {
    long seed = std::stol(argv[1]);
    double side = std::stod(argv[2]);
    long ncside = std::stol(argv[3]);
    long long n_part = std::stoll(argv[4]);
    long long time_steps = std::stoll(argv[5]);

    if (rank > ncside - 1) {
      MPI_Finalize();
      return 0;
    }

    init_particles(seed, side, ncside, n_part, particles, rank);

    exec_time = -omp_get_wtime();
    PartialGrid grid = init_grid(rank, nprocs, side, ncside, particles);
    for (long long ll = 0; ll < time_steps; ll++)
      grid.update();
    grid.sync_final_state();
    exec_time += omp_get_wtime();

    if (rank == 0) {
      fprintf(stderr, "%.1fs\n", exec_time);
      print_result(grid); // to stdout
    }
  } catch (const std::exception &e) {
    std::cerr << "Rank: " << rank << "\n";
    std::cerr << "Error: " << e.what() << "\n";
    MPI_Finalize();
    return 1;
  }

  MPI_Finalize();
  return 0;
}
