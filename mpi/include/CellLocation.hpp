#pragma once

class CellLocation {
public:
  long rank;
  double x, y;
  bool counted;
  CellLocation(double x, double y, long rank)
      : rank(rank), x(x), y(y), counted(false) {}
};
