#ifndef __H2D_TEST_CAND_PROJ_H
#define __H2D_TEST_CAND_PROJ_H

//test case
class TestCase {
  int _func_quad_order; ///< Function order.
  int _start_quad_order; ///< Start order of an element.
  double** _poly_matrix; ///< MxN matrix of polygonal coefficients. The value is calculates [x^0, ..., x^(M-1)] * matrix * [y^0, ..., y^(N-1)]^T
public:
  TestCase(int func_quad_order);
  ~TestCase();

  bool should_match(const RefinementSelectors::OptimumSelector::Cand& cand); ///< Returns true if the refinement should match the function.

  std::string title() const;
  int quad_order() const { return _func_quad_order; }
  int start_quad_order() const { return _start_quad_order; }
  double** poly_matrix() { return _poly_matrix; };
};

//text functions
extern scalar func_val(double** poly_matrix, int quad_order, double x, double y);
extern scalar func_dx(double** poly_matrix, int quad_order, double x, double y);
extern scalar func_dy(double** poly_matrix, int quad_order, double x, double y);

#endif
