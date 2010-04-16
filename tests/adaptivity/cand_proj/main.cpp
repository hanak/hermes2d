#include <hermes2d.h>
#include <solver_umfpack.h>
#include "functions.h"

using namespace RefinementSelectors;

// This test tests projection of a candidate in H1 space on a quad.
/* global definitions */
#undef ERROR_SUCCESS
#undef ERROR_FAILURE
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
#define H2D_TEST_P_INIT 1 /* initial order on a mesh */
#define H2D_TEST_ELEM_ID 0 /* ID of an alement which is going to be handled */
#define H2D_TEST_ZERO 1e-12 /* Numerical zero for double datatype */
#define delete_not_null(__ptr) if (__ptr != NULL) delete __ptr;

enum FailCode { ///< Failure code: used to show failure of tests.
  H2D_TEST_NOT_DONE,
  H2D_TEST_FAILED,
  H2D_TEST_SUCCESS
};

/* global variables */
Mesh* mesh = NULL;
H1Shapeset* shapeset = NULL;
PrecalcShapeset* pss = NULL;
H1Space* space = NULL;
WeakForm* weakform = NULL;
UmfpackSolver* solver = NULL;
TestCase* cur_test_case = NULL;

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
    space->set_uniform_order(H2D_TEST_P_INIT);
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
void show_fail_matrix(FailCode** fail_matrix, const int max_quad_order) {
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
      char fail_char = '?';
      switch (fail_matrix[i][k]) {
        case H2D_TEST_NOT_DONE: fail_char = '-'; break;
        case H2D_TEST_FAILED: fail_char = 'X'; break;
        case H2D_TEST_SUCCESS: fail_char = ' '; break;
      }
      for(int j = 0; j < NUMBER_W; j++)
        str << fail_char;
    }

    //print row
    info(" %s", str.str().c_str());
  }
}

/// Test
int test(bool tri) {
   bool failed = false;

  //prepare selector
  H1ProjBasedSelector selector(H2D_HP_ANISO, 1.0, H2DRS_DEFAULT_ORDER, shapeset);

  //prepare cases
  OrderPermutator order_perm(H2D_MAKE_QUAD_ORDER(0, 0), H2D_MAKE_QUAD_ORDER(H2DRS_MAX_ORDER, H2DRS_MAX_ORDER), tri);
  const int min_quad_order = H2D_MAKE_QUAD_ORDER(2, 2);
  //OrderPermutator order_perm(H2D_MAKE_QUAD_ORDER(2, 3), H2D_MAKE_QUAD_ORDER(2, 3), tri);

  //prepare place for result summary
  FailCode** fail_matrix = new_matrix<FailCode>(H2D_GET_V_ORDER(order_perm.get_end_quad_order())+1, H2D_GET_H_ORDER(order_perm.get_end_quad_order())+1);

  //process cases
  do {
    FailCode test_result;
    if (order_perm.get_order_h() >= H2D_GET_H_ORDER(min_quad_order) && order_perm.get_order_v() >= H2D_GET_V_ORDER(min_quad_order)) {
      test_result = H2D_TEST_SUCCESS;

      //prepare test case
      TestCase test_case(order_perm.get_quad_order());
      cur_test_case = &test_case;
      verbose("!test case: %s", test_case.title().c_str());

      //process
      space->set_element_order(H2D_TEST_ELEM_ID, test_case.start_quad_order());
      int ndofs = space->assign_dofs();

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

      //select candidate
      ElementToRefine refinement;
      Element* e = mesh->get_element(H2D_TEST_ELEM_ID);
      int order = space->get_element_order(H2D_TEST_ELEM_ID);
      selector.select_refinement(e, order, &rsln, refinement);

      //check candidates
      const std::vector<OptimumSelector::Cand>& candidates = selector.get_candidates();
      std::vector<OptimumSelector::Cand>::const_iterator cand = candidates.begin();
      while (cand != candidates.end()) {
        if (cur_test_case->should_match(*cand) && abs(cand->error) > H2D_TEST_ZERO) {
          std::stringstream str;
          str << *cand;
          verbose(" invalid candidate: %s", str.str().c_str());
          test_result = H2D_TEST_FAILED;
        }
        cand++;
      }
      if (test_result == H2D_TEST_SUCCESS)
        verbose("test success");
      else {
        verbose("test failed!");
        failed = true;
      }
    }
    else
      test_result = H2D_TEST_NOT_DONE;
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
