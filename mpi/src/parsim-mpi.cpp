#include "Cell.hpp"
#include "CellInfo.hpp"
#include "PartialCell.hpp"
#include "PartialGrid.hpp"
#include "Particle.hpp"
#include <unordered_map>
#include <unordered_set>
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
    par[i].x = rnd01() * side;
    par[i].y = rnd01() * side;
    par[i].vx = (rnd01() - 0.5) * side / ncside / 5.0;
    par[i].vy = (rnd01() - 0.5) * side / ncside / 5.0;

    par[i].m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
  }
  par[0].first_particle = true;
}

void partition_grid(int nprocs, long ncside, long n_part, double cell_size,
                    std::vector<CellInfo> &partition) {
  partition.reserve(ncside * ncside);
  // Process we are currently assigning cells to
  int proc = 0;
  // Weight each proc is responsible for
  std::vector<long> weight(nprocs, 0);
  // We want each proc to have roughly the same weight
  long target = n_part / nprocs;
  int size = std::max(5, (int)std::sqrt((ncside * ncside) / nprocs));
  std::pair<long, long> loc = {0, 0};

  while (loc.first != ncside && loc.second != ncside) {
    long xsize = loc.first;
    long ysize = std::min(loc.second + size, ncside);

    for (; weight[proc] < target && xsize < ncside; ++xsize) {
      for (long i = loc.second; i < ysize; ++i) {
        partition[xsize + i * ncside].rank = proc;
        weight[proc] += partition[xsize + i * ncside].weight;
      }
    }

    if (xsize == ncside && ysize != ncside) {
      loc.first = 0;
      loc.second = ysize;
    } else {
      loc.first = xsize;
    }
    proc = (proc + 1) % nprocs;
  }
}

PartialGrid init_grid(int rank, int nprocs, double side, long ncside,
                      std::vector<Particle> &pv) {
  double cell_size = side / ncside;
  std::vector<CellInfo> partition;

  // Creating cell information
  for (long i = 0; i < ncside; ++i) {
    double y = i * cell_size;
    for (long j = 0; j < ncside; ++j) {
      double x = j * cell_size;
      partition.emplace_back(j + i * ncside, x, y);
    }
  }

  // Assign particles to corresponding cells
  for (auto &p : pv) {
    long cell_x = static_cast<long>(p.x / cell_size);
    long cell_y = static_cast<long>(p.y / cell_size);
    long idx = cell_x + cell_y * ncside;
    partition[idx].add_particle(p);
  }

  partition_grid(nprocs, ncside, pv.size(), cell_size, partition);
  PartialGrid grid(rank, nprocs, side, ncside);

  for (auto &info : partition) {
    int cid = info.id;
    std::unordered_set<int> adj;
    bool owner_is_adjacent = false;
    std::vector<long> adj_cells = PartialGrid::get_adjacent_cells(cid, ncside);

    for (const auto &id : adj_cells) {
      if (partition[id].rank == rank) {
        owner_is_adjacent = true;
        continue;
      }

      if (info.rank != rank) {
        if (owner_is_adjacent) {
          break;
        } else {
          continue;
        }
      }

      adj.insert(partition[id].rank);
      grid.add_adjacent_rank(partition[id].rank);
      if (!partition[id].counted) {
        partition[id].counted = true;
        grid.increment_adjacent_rank(partition[id].rank);
      }
    }

    if (info.rank == rank) {
      Cell cell = Cell(info, adj, cell_size);
      grid.add_local_cell(cell);
    } else if (owner_is_adjacent) {
      PartialCell cell = PartialCell(info, cell_size);
      grid.add_adjacent_cell(cell);
    }
  }

  // Initialize cells
  for (auto &it : grid.fully_local_cells)
    it.second.check_collisions();
  for (auto &it : grid.partially_local_cells)
    it.second.check_collisions();

  return grid;
}

void print_result(PartialGrid &g) {
  Particle pf = g.get_first_particle();
  fprintf(stdout, "%.3f %.3f\n", pf.x, pf.y);

  std::cout << g.get_collisions() << std::endl;
}

void simulation(PartialGrid &grid, long long time_steps) {
  grid.update(time_steps);
  grid.sync_final_state();
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

    init_particles(seed, side, ncside, n_part, particles);

    exec_time = -omp_get_wtime();
    PartialGrid grid = init_grid(rank, nprocs, side, ncside, particles);
    simulation(grid, time_steps);
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
