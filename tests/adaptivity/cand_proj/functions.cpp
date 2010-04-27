#include <hermes2d.h>
#include "functions.h"

using namespace RefinementSelectors;

TestCase::TestCase(int func_quad_order)
  : m_func_quad_order(func_quad_order)
  , m_start_quad_order(make_quad_order(std::max(0, get_h_order(func_quad_order)-1), std::max(0, get_v_order(func_quad_order)-1))) {
  //create and fill poly matrix
  const int order_h = get_h_order(func_quad_order), order_v = get_v_order(func_quad_order);
  m_poly_matrix = new_matrix<double>(order_v + 1, order_h + 1);
  for(int i = 0; i <= order_v; i++) {
    for(int k = 0; k <= order_h; k++) {
      int rnd_val = (rand() - RAND_MAX/2);
      if (rnd_val == 0)
        rnd_val = 1;
      m_poly_matrix[i][k] = (rnd_val * 2.0) / RAND_MAX;
      ////DEBUG-BEGIN
      //if (k == 0)
      //  m_poly_matrix[i][k] = 1;
      //else
      //  m_poly_matrix[i][k] = 0;
      ////DEBUG-END
    }
  }
}

TestCase::~TestCase() { delete[] m_poly_matrix; };

bool TestCase::should_match(const RefinementSelectors::OptimumSelector::Cand& cand) {
  int order_h = get_h_order(m_func_quad_order), order_v = get_v_order(m_func_quad_order);
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
  str << "Polynom order H:" << get_h_order(m_func_quad_order) << "; V:" << get_v_order(m_func_quad_order);
  return str.str();
}

scalar func_val(double** poly_matrix, int quad_order, double x, double y) {
  int order_h = get_h_order(quad_order), order_v = get_v_order(quad_order);
  double res = 0;

  //x * M * y
  for(int i = 0; i <= order_v; i++) {
    double val = 0;
    for(int k = 0; k <= order_h; k++) {
      val += pow(x, k) * poly_matrix[i][k];
    }
    res += val * pow(y, i);
  }

  return res;
}

scalar func_dx(double** poly_matrix, int quad_order, double x, double y) {
  int order_h = get_h_order(quad_order), order_v = get_v_order(quad_order);
  double res = 0;

  //dx/dx * M * y
  for(int i = 0; i <= order_v; i++) {
    double val = 0;
    for(int k = 1; k <= order_h; k++) {
      val += pow(x, k-1) * k * poly_matrix[i][k];
    }
    res += val * pow(y, i);
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
    for(int k = 0; k <= order_h; k++) {
      val += pow(x, k) * poly_matrix[i][k];
    }
    res += val * i * pow(y, i-1);
  }

  return res;
}
