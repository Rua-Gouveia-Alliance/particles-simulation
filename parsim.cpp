#define _USE_MATH_DEFINES
#include <math.h>
#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1
#include <vector>

class Particle {
  double _x;
  double _y;
  double _vx;
  double _vy;
  double _m;

public:
  Particle(double x, double y, double vx, double vy, double m)
      : _x(x), _y(y), _vx(vx), _vy(vy), _m(m) {}
};

unsigned int seed;

void init_r4uni(int input_seed) { seed = input_seed + 987654321; }

double rnd_uniform01() {
  int seed_in = seed;
  seed ^= (seed << 13);
  seed ^= (seed >> 17);
  seed ^= (seed << 5);
  return 0.5 + 0.2328306e-09 * (seed_in + (int)seed);
}

double rnd_normal01() {
  double u1, u2, z, result;
  do {
    u1 = rnd_uniform01();
    u2 = rnd_uniform01();
    z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    result = 0.5 + 0.15 * z; // Shift mean to 0.5 and scale
  } while (result < 0 || result >= 1);
  return result;
}

// TODO: Troquei um bocado esta funcao, ns se podemos
void init_particles(long seed, double side, long ncside, long long n_part,
                    std::vector<Particle> &par /* particle par */) {
  double (*rnd01)() = rnd_uniform01;
  long long i;
  if (seed < 0) {
    rnd01 = rnd_normal01;
    seed = -seed;
  }
  init_r4uni(seed);
  for (i = 0; i < n_part; i++) {
    double x = rnd01() * side;
    double y = rnd01() * side;
    double vx = (rnd01() - 0.5) * side / ncside / 5.0;
    double vy = (rnd01() - 0.5) * side / ncside / 5.0;
    double m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
    par.emplace_back(x, y, vx, vy, m);
    // par[i].x = rnd01() * side;
    // par[i].y = rnd01() * side;
    // par[i].vx = (rnd01() - 0.5) * side / ncside / 5.0;
    // par[i].vy = (rnd01() - 0.5) * side / ncside / 5.0;
    // par[i].m = rnd01() * 0.01 * (ncside * ncside) / n_part / G * EPSILON2;
  }
}
