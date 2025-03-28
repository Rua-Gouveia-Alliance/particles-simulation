#pragma once

class CellLocation {
public:
  double x, y;
  long rank;
  CellLocation(double x, double y, long rank) : rank(rank), x(x), y(y) {}
};
