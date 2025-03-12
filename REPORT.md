# Project Report

Group 5:
- Gonçalo Rua 102604
- Rodrigo Arêde 102606
- João Gouveia 102611

## 1. Introduction

This stage of the project focuses on parallelizing the serial version of a particles simulation, using OpenMP.
This simulation models the gravitational interaction between particles in a 2D space, and the goal is to improve
the computational efficiency through parallelization.

## 2. Approach

The 2D space is divided into a grid, where each cell's center of mass is used to approximate the gravitational interactions.
We took the approach of parallelizing the cells of the grid, with the objective of each cell being processed concurrently with
other cells. We also parallelized the generation of the initial grid which required calculations to attribute each particle to
the correspoding cell where it's initially located.

## 3. Synchronization concerns

We had some challenges with the synchronization of the different threads, specifically when we encountered problems with data
dependencies between them. One clear example was the appending of newly updated particles to the cells, while other particles from
different cells were being computed. To solve this problem we had to separate this stage in two, which led us to temporarly save
all the new particles in a `vector<vector<Particles>>`, and only after all the particles in all the cells are done computing we
proceed to add them to the correct cells.

## 4. Load Balancing

#TODO With OpenMP we used dynamic scheduling to evenly distrubte computation across threads. This technique distributes workload among threads at runtime rather than assigning it statically before execution. Comparing the static vs dynamic approach, the latter gives us the best performance of the two.

## 5. Performance Results

For analysing the performance of our implementation we used the following inputs:

| Test Number | Input                  |
|-------------|------------------------|
| 1           | 12672 0.05 3 10 10     |
| 2           | 5893 0.05 3 10 10      |
| 3           | 8555 0.05 3 10 10      |
| 4           | 12 100 5 10000 10000   |
| 5           | -11 3500 20 500000 10  |
| 6           | 1 5000 100 1000000 4   |
| 7           | 1 5000 100 1000000 100 |
| 8           | 1 5000 20 1000000 10   |
| 9           | 1 1000 3 10000 10000   |
| 10          | 3 5000 50 1000000 300  |
| 11          | 3 5000 50 1000000 500  |
| 12          | -1 1000 30 100000 1000 |

### 5.1 Table

The table below shows the execution times (in seconds) for the serial and OpenMP versions of the program:

| Test Number | 1 thread (serial) | 2 threads | 4 threads | 8 threads |
|-------------|-------------------|-----------|-----------|-----------|
| 1           | 0.0s              | 0.0s      | 0.0s      | 0.0s      |
| 2           | 0.0s              | 0.0s      | 0.0s      | 0.0s      |
| 3           | 0.0s              | 0.0s      | 0.0s      | 0.0s      |
| 4           | 53.9s             | 32.1s     | 21.8s     | 19.1s     |
| 5           | 49.7s             | 26.4s     | 25.2s     | 20.0s     |
| 6           | 1.6s              | 1.0s      | 0.8s      | 0.8s      |
| 7           | 34.9s             | 22.7s     | 17.5s     | 17.0s     |
| 8           | 57.1s             | 29.1s     | 15.6s     | 9.4s      |
| 9           | 245.9s            | 145.0s    | 93.0s     | 74.3s     |
| 10          | 289.4s            | 163.4s    | 100.9s    | 81.3s     |
| 11          | 481.9s            | 268.1s    | 167.9s    | 133.5s    |
| 12          | 81.8s             | 46.9s     | 46.2s     | 40.2s     |

## 6. Conclusion

The project successfully implemented a parallel particles simulation, demonstrating significant speedups using OpenMP. The next step is to use MPI to parallelize the workload across different machines, and exploring an hybrid MPI+OpenMP approach for even better scalability.
