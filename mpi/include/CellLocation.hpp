#pragma once

class CellLocation {
public:
  long rank;
  double x, y;
  bool counted;
  CellLocation(double x, double y, long rank)
      : x(x), y(y), rank(rank), counted(false) {}
};
