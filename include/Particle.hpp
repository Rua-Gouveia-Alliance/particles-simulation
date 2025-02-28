#ifndef __PARTICLE_HPP__
#define __PARTICLE_HPP__

class Particle {
private:
public:
  double x, y, vx, vy, m;
  bool first_particle;

  Particle();
  Particle(double m);
  Particle(double x, double y, double vx, double vy, double m);

  void print_info() const;
};

#endif
