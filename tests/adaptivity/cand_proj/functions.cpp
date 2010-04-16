#include <hermes2d.h>
#include "functions.h"

using namespace RefinementSelectors;

TestCase::TestCase(int func_quad_order)
  : _func_quad_order(func_quad_order), _start_quad_order(make_quad_order(get_h_order(func_quad_order)-1, get_v_order(func_quad_order)-1)) {
  //create and fill poly matrix
  const int order_h = get_h_order(func_quad_order), order_v = get_v_order(func_quad_order);
  _poly_matrix = new_matrix<double>(order_v + 1, order_h + 1);
  for(int i = 0; i <= order_v; i++) {
    for(int k = 0; k <= order_h; k++) {
      int rnd_val = (rand() - RAND_MAX/2);
      if (rnd_val == 0)
        rnd_val = 1;
      if (i != order_v || k != order_h)
        rnd_val = 0;
      _poly_matrix[i][k] = (rnd_val * 2.0) / RAND_MAX;
    }
  }
}

TestCase::~TestCase() { delete[] _poly_matrix; };

bool TestCase::should_match(const OptimumSelector::Cand& cand) { ///< Returns true if the refinement should match the function.
  int order_h = get_h_order(_func_quad_order), order_v = get_v_order(_func_quad_order);
  int num_sons = cand.get_num_sons();
  for(int i = 0; i < num_sons; i++) {
    int son_order_h = get_h_order(cand.p[i]), son_order_v = get_v_order(cand.p[i]);
    if (son_order_h < order_h || son_order_v < order_v)
      return false;
  }
  return true;
}

std::string TestCase::title() const {
  std::stringstream str;
  str << "Polynom order H:" << get_h_order(_func_quad_order) << "; V:" << get_v_order(_func_quad_order);
  return str.str();
}

scalar func_val(double** poly_matrix, int quad_order, double x, double y) {
  int order_h = get_h_order(quad_order), order_v = get_v_order(quad_order);
  double res = 0;

  //x * M * y
  double y_pow = 1;
  for(int i = 0; i <= order_v; i++) {
    double val = 0;
    double x_pow = 1;
    for(int k = 0; k <= order_h; k++) {
      val += x_pow * poly_matrix[i][k];
      x_pow *= x;
    }
    res += val * y_pow;
    y_pow *= y;
  }

  return res;
}

scalar func_dx(double** poly_matrix, int quad_order, double x, double y) {
  int order_h = get_h_order(quad_order), order_v = get_v_order(quad_order);
  double res = 0;

  //dx/dx * M * y
  double y_pow = 1;
  for(int i = 0; i <= order_v; i++) {
    double val = 0;
    double x_pow = 1;
    for(int k = 1; k <= order_h; k++) {
      val += x_pow * k * poly_matrix[i][k];
      x_pow *= x;
    }
    res += val * y_pow;
    y_pow *= y;
  }

  return res;
}

scalar func_dy(double** poly_matrix, int quad_order, double x, double y) {
  int order_h = get_h_order(quad_order), order_v = get_v_order(quad_order);
  double res = 0;

  //x * M * dy/dy
  double y_pow = 1;
  for(int i = 1; i <= order_v; i++) {
    double val = 0;
    double x_pow = 1;
    for(int k = 0; k <= order_h; k++) {
      val += x_pow * poly_matrix[i][k];
      x_pow *= x;
    }
    res += val * y_pow * i;
    y_pow *= y;
  }

  return res;
}
