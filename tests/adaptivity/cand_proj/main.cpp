/// \addtogroup t_cand_proj Candidate projection
/// \{
/// \brief This test tests projection of a candidate in H1 space on a quad.
///
/// For each combination of orders it creates a polynom \f$ p(x,y) = \sum_h \sum_v a_{h,v} x^h y^v \f$.
/// The test ensures that all coefficients \f$a_{h,v}\f$ are not zero. The polynom is defined
/// in the reference domain.
/// Then, it creates a mesh which consist of a single element. The order of the element is
/// \f$(h-1, v-1)\f$ where \f$h\f$ is a current horizontal order and \f$v\f$ is current vertical order.
/// Using the mesh and the polynom it calculates a solution and reference solution. The reference solution
/// is used the control the selection of candidates. Using projection based selector,
/// the test generates HP_ANISO candidates and calculates their errors.
///
/// The test succeeds if:
///  - Values of reference solution equals to values of the original function in random points.
///  - Errors of all candidates whose orders of all sons are greater or equal to order of the polynom are zero.
///
/// Output table legend:
///  - '-': Test was not done
///  - ' ': Test succesfull
///  - 'S': Values of the reference solution differs from values of the test function despite that the order of the reference solution is either the same of higher than the test function.
///  - 'C': Some candidates has non-zero error despite that order of their sons is greater or equal to the order of the test function.
///
/// \dontinclude hermes2d.h
/// \dontinclude solver_umfpack.h

#include <hermes2d.h>
#include <solver_umfpack.h>
#include "functions.h"

using namespace RefinementSelectors;

/* global definitions */
#undef ERROR_SUCCESS
#undef ERROR_FAILURE
#define ERROR_SUCCESS 0 ///< Test return code if success.
#define ERROR_FAILURE -1 ///< Test return code if fails.
#define H2D_TEST_ELEM_ID 0 ///< ID of an alement which is going to be handled.
#define H2D_TEST_ZERO 2e-12 ///< Numerical zero. Since polynoms are defined on a reference domain, some high-order polynoms might yield error a little bit higher than \f$10^{-12}\f$.
#define delete_not_null(__ptr) if (__ptr != NULL) delete __ptr; ///< Deletes an instance if the pointer is not NULL.

#define H2D_TEST_NOT_DONE -1 ///< Flag: Combination of orders was not tested.
#define H2D_TEST_SUCCESS 0x00 ///< Flag: Combination of orders was tested successfully.
#define H2D_CAND_FAILED 0x01  ///< Flag: Appropriate candidates have non-zero error.
#define H2D_RSLN_FAILED 0x02  ///< Flag: failure of the condition of equality between ref. solution failed and the polynom.

/* global variables */
Mesh* mesh = NULL; ///< Mesh used by the test.
H1Shapeset* shapeset = NULL; ///< Shapeset used by the test.
PrecalcShapeset* pss = NULL; ///< Precalculated shapeset used by the test.
H1Space* space = NULL; ///< Space used by the test.
WeakForm* weakform = NULL; ///< Weakform used by the test.
UmfpackSolver* solver = NULL; ///< Solved used by the test.
TestCase* cur_test_case = NULL; ///< Current test case: required by callbacks.

/// Boundary condition types
int bc_types(int marker) { return BC_NONE; }

/// Bilinear form.
template<typename T, typename D>
D h1_biform(int n, double *wt, Func<T> *u, Func<T> *v, Geom<T> *e, ExtData<D> *data) {
  return int_u_v<T, T>(n, wt, u, v) + int_grad_u_grad_v<T, T>(n, wt, u, v);
}

/// Linear form: order
Ord h1_liform(int point_cnt, double *weights, Func<Ord> *values_v, Geom<Ord> *geometry, ExtData<Ord> *values_fnc_ext) {
  return Ord(H2D_GET_H_ORDER(cur_test_case->quad_order()) + H2D_GET_V_ORDER(cur_test_case->quad_order()) + 2*values_v->val->get_order());
}

/// Linear form: value
scalar h1_liform(int point_cnt, double *weights, Func<double> *values_v, Geom<double> *geometry, ExtData<scalar> *values_fnc_ext) {
  scalar result = 0;
  for(int i = 0; i < point_cnt; i++) {
    double x = geometry->x[i], y = geometry->y[i];
    result += weights[i] * (values_v->val[i] * func_val(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y)
      + values_v->dx[i] * func_dx(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y)
      + values_v->dy[i] * func_dy(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y));
  }
  return result;
}

/// Function that is used to exact solution.
scalar exact_func(double x, double y, scalar& dx, scalar& dy) {
  scalar value = func_val(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y);
  dx = func_dx(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y);
  dy = func_dy(cur_test_case->poly_matrix(), cur_test_case->quad_order(), x, y);
  return value;
}

/// Initialize internal data structures.
bool init(bool tri) {
  try {
    // shapeset and cache
    shapeset = new H1Shapeset();
    pss = new PrecalcShapeset(shapeset);

    // mesh (1 element)
    mesh = new Mesh();
    if (tri) {
      const int vertex_num = 3, tria_num = 1, quad_num = 0, marker_num = 3;
      double2 vertex_array[3] = {{-1,-1}, { 1,-1}, {-1, 1}};
      int4 tria_array[1] = {{0, 1, 2, 0}};
      int5 *quad_array = NULL;
      int3 marker_array[3] = {{0, 1, 2}, {1, 2, 1}, {2, 0, 3}};
      mesh->create(vertex_num, vertex_array, tria_num, tria_array, quad_num, quad_array, marker_num, marker_array);
    }
    else {
      const int vertex_num = 4, tria_num = 0, quad_num = 1, marker_num = 4;
      double2 vertex_array[4] = {{-1,-1}, { 1,-1}, { 1, 1}, {-1, 1}};
      int4 *tria_array = NULL;
      int5 quad_array[1] = {{0, 1, 2, 3, 0}};
      int3 marker_array[4] = {{0, 1, 3}, {1, 2, 2}, {2, 3, 4}, {3, 0, 1}};
      mesh->create(vertex_num, vertex_array, tria_num, tria_array, quad_num, quad_array, marker_num, marker_array);
    }

    // finite element space
    space = new H1Space(mesh, shapeset);
    space->set_bc_types(bc_types);
    space->set_uniform_order(1);
    space->assign_dofs();

    // weakform
    weakform = new WeakForm(1);
    weakform->add_biform(0, 0, callback(h1_biform), H2D_SYM);
    weakform->add_liform(0, h1_liform, h1_liform, H2D_ANY);

    //solver
    solver = new UmfpackSolver();

    return true;
  }
  catch (std::exception& e) {
    info("failed: %s", e.what());
    return false;
  }
}

/// Cleanup
void cleanup() {
  delete_not_null(solver);
  delete_not_null(weakform);
  delete_not_null(space);
  delete_not_null(mesh);
  delete_not_null(pss);
  delete_not_null(shapeset);
}

/// Prints failure matrix
void show_fail_matrix(int** fail_matrix, const int max_quad_order) {
  const int max_order_h = H2D_GET_H_ORDER(max_quad_order), max_order_v = H2D_GET_V_ORDER(max_quad_order);
#define NUMBER_W 2
#define NUMBER_FMT "%02d"
  char buff_number[1024];

  info("!Test summary (V/H)");

  //print header
  {
    std::stringstream str;
    for(int i = 0; i < NUMBER_W; i++)
      str << ' ';
    for(int i = 0; i <= max_order_h; i++) {
      str << '|';
      sprintf(buff_number, NUMBER_FMT, i);
      str << buff_number;
    }
    info(" %s", str.str().c_str());
  }

  //print body
  for(int i = 0; i <= max_order_v; i++) {
    //build row head
    std::stringstream str;
    sprintf(buff_number, NUMBER_FMT, i);
    str << buff_number;

    //build row body
    for(int k = 0; k <= max_order_h; k++) {
      str << '|';
      if (fail_matrix[i][k] == H2D_TEST_NOT_DONE) {
        for(int j = 0; j < NUMBER_W; j++)
          str << '-';
      }
      else {
        if ((fail_matrix[i][k] & H2D_RSLN_FAILED) != 0)
          str << 'S';
        else
          str << ' ';
        if ((fail_matrix[i][k] & H2D_CAND_FAILED) != 0)
          str << 'C';
        else
          str << ' ';
        for(int j = 2; j < NUMBER_W; j++)
          str << ' ';
      }
    }

    //print row
    info(" %s", str.str().c_str());
  }
}

/// Test reference solution whether its values matche values of the polynom in the test case.
bool test_ref_solution(Solution& rsln) {
  verbose("Testing reference solution");

  //check FN order
  int test_quad_order = cur_test_case->quad_order();

  //check function values
  bool failed = false;
  Mesh* mesh = rsln.get_mesh();
  Element* element;
  for_all_active_elements(element, mesh) {
    //set element
    rsln.set_active_element(element);
    int rsln_order = rsln.get_fn_order();
    if (rsln_order < H2D_GET_H_ORDER(test_quad_order) || rsln_order < H2D_GET_V_ORDER(test_quad_order)) {
      verbose(" Failed: reference solution order (%d/%d) is lower than test case order (%d/%d) at element #%d", rsln_order, rsln_order, H2D_GET_H_ORDER(test_quad_order), H2D_GET_V_ORDER(test_quad_order), element->id);
      return true;
    }

    //set GIP order at which points will be inspected
    int check_gip_order = 2*rsln_order;
    rsln.set_quad_order(check_gip_order);

    //get rsln values
    scalar* rsln_vals = rsln.get_fn_values(0);

    //get physical locations
    RefMap* ref_map = rsln.get_refmap();
    scalar* phys_x = ref_map->get_phys_x(check_gip_order);
    scalar* phys_y = ref_map->get_phys_y(check_gip_order);

    //compare to test-case values
    const int num_pt = rsln.get_quad_2d()->get_num_points(check_gip_order);
    for(int i = 0; i < num_pt; i++) {
      double func_value = func_val(cur_test_case->poly_matrix(), test_quad_order, phys_x[i], phys_y[i]);
      if (abs(func_value - rsln_vals[i]) > H2D_TEST_ZERO) {
        verbose(" Failed: rsln and function value differs");
        return true;
      }
    }
  }

  return false;
}

/// Test
int test(bool tri) {
   bool failed = false;

  //prepare selector
  H1ProjBasedSelector selector(H2D_HP_ANISO, 1.0, H2DRS_DEFAULT_ORDER, shapeset);

  //prepare cases
  int min_quad_order = -1, max_quad_order = -1;
  if (tri) {
    min_quad_order = H2D_MAKE_QUAD_ORDER(2, 0);
    max_quad_order = H2D_MAKE_QUAD_ORDER(H2DRS_MAX_ORDER, 0);
    //min_quad_order = H2D_MAKE_QUAD_ORDER(4, 0); //DEBUG
    //max_quad_order = H2D_MAKE_QUAD_ORDER(5, 0); //DEBUG
  }
  else {
    min_quad_order = H2D_MAKE_QUAD_ORDER(2, 2);
    max_quad_order = H2D_MAKE_QUAD_ORDER(H2DRS_MAX_ORDER, H2DRS_MAX_ORDER);
  }
  OrderPermutator order_perm(min_quad_order, max_quad_order, false);

  //prepare place for result summary
  int end_order_v = H2D_GET_V_ORDER(order_perm.get_end_quad_order());
  int end_order_h = H2D_GET_H_ORDER(order_perm.get_end_quad_order());
  int** fail_matrix = new_matrix<int>(end_order_v+1, end_order_h+1);
  for(int i = 0; i <= end_order_v; i++)
    for(int k = 0; k <= end_order_h; k++)
      fail_matrix[i][k] = H2D_TEST_NOT_DONE;

  //process cases
  do {
    int test_result = H2D_TEST_SUCCESS;

    //prepare test case
    TestCase test_case(order_perm.get_quad_order());
    cur_test_case = &test_case;
    verbose("!test case: %s", test_case.title().c_str());

    //process
    space->set_element_order(H2D_TEST_ELEM_ID, test_case.start_quad_order());
    int ndofs = space->assign_dofs();
    space->assign_dofs();

    //create and solve the reference system
    Solution sln, rsln;
    LinSystem ls(weakform, solver);
    ls.set_spaces(1, space);
    ls.set_pss(1, pss);
    ls.assemble();
    ls.solve(1, &sln);
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &rsln);

    ////DEBUG-BEGIN: display
    //ScalarView view;
    //ExactSolution esln(mesh, exact_func);
    //DiffFilter diff_filter(&sln, &esln);
    //view.show(&diff_filter);
    //View::wait(H2DV_WAIT_KEYPRESS);
    ////DEBUG-END

    //check projected functions
    if (test_ref_solution(rsln))
      test_result |= H2D_RSLN_FAILED;

    //select candidate
    ElementToRefine refinement;
    Element* e = mesh->get_element(H2D_TEST_ELEM_ID);
    int order = space->get_element_order(H2D_TEST_ELEM_ID);
    selector.select_refinement(e, order, &rsln, refinement);

    //check candidates
    verbose("Testing candidates");
    const std::vector<OptimumSelector::Cand>& candidates = selector.get_candidates();
    std::vector<OptimumSelector::Cand>::const_iterator cand = candidates.begin();
    while (cand != candidates.end()) {
      if (cur_test_case->should_match(*cand) && abs(cand->error) > H2D_TEST_ZERO) {
        std::stringstream str;
        str << *cand;
        verbose(" Invalid candidate: %s", str.str().c_str());
        test_result |= H2D_CAND_FAILED;
      }
      cand++;
    }
    if (test_result == H2D_TEST_SUCCESS)
      verbose("Test success");
    else {
      verbose("Test failed!");
      failed = true;
    }
    fail_matrix[order_perm.get_order_v()][order_perm.get_order_h()] = test_result;
  } while(order_perm.next());

  //print result matrix
  show_fail_matrix(fail_matrix, order_perm.get_end_quad_order());

  //clenup
  delete[] fail_matrix;

  if (failed)
    return ERROR_FAILURE;
  else
    return ERROR_SUCCESS;
}

/// Test entry-point.
int main(int argc, char* argv[]) {
  //check if triangle is requested
  bool tri = false;
  if (argc > 1 && strcmp(argv[1], "-tri") == 0)
    tri = true;

  //initialize
  int result = ERROR_FAILURE;
  if (init(tri)) {
    result = test(tri);
  }
  cleanup();

  if (result == ERROR_SUCCESS)
  {
    info("!Test: Success");
  }
  else {
    info("!Test: Failed!");
  }
  return result;
}

/// \}
