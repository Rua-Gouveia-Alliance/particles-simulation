#pragma once

class Mass {
private:
  long _id;
  double _x, _y, _val;

public:
  Mass(long id) : _id(id), _x(-1), _y(-1), _val(0) {}

  Mass(long id, double x, double y, double val)
      : _id(id), _x(x), _y(y), _val(val) {}

  long id() const { return _id; }
  double x() const { return _x; }
  double y() const { return _y; }
  double val() const { return _val; }

  void x(double x) { _x = x; }
  void y(double y) { _y = y; }
  void val(double val) { _val = val; }
};
