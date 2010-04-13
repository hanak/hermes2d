#include "../common.h"
#include "../matrix.h"
#include "../solution.h"
#include "../shapeset_h1_all.h"
#include "../element_to_refine.h"
#include "h1_proj_based_selector.h"

namespace RefinementSelectors {
  H1Shapeset H1ProjBasedSelector::default_shapeset;

  H1ProjBasedSelector::H1ProjBasedSelector(CandList cand_list, double conv_exp, int max_order, H1Shapeset* user_shapeset)
    : ProjBasedSelector(cand_list, conv_exp, max_order, user_shapeset == NULL ? &default_shapeset : user_shapeset, Range<int>(1,1), Range<int>(2, H2DRS_MAX_ORDER)) {}

  void H1ProjBasedSelector::set_current_order_range(Element* element) {
    current_max_order = this->max_order;
    if (current_max_order == H2DRS_DEFAULT_ORDER)
      current_max_order = (20 - element->iro_cache)/2 - 1; // default
    else
      current_max_order = std::min(current_max_order, (20 - element->iro_cache)/2 - 1); // user specified
    current_min_order = 1;
  }

  scalar** H1ProjBasedSelector::precalc_ref_solution(int inx_son, Solution* rsln, Element* element, int intr_gip_order) {
    //set element and integration order
    rsln->set_active_element(element);
    rsln->set_quad_order(intr_gip_order);

    //fill with values
    scalar** rvals_son = precalc_rvals[inx_son];
    rvals_son[H2D_FN_VALUE] = rsln->get_fn_values(0);
    rvals_son[H2D_FN_DX] = rsln->get_dx_values(0);
    rvals_son[H2D_FN_DY] = rsln->get_dy_values(0);

    return rvals_son;
  }

  double** H1ProjBasedSelector::build_projection_matrix(Shapeset& shapeset,
    double3* gip_points, int num_gip_points,
    const int* shape_inx, const int num_shapes) {
    //allocate
    double** matrix = new_matrix<double>(num_shapes, num_shapes);

    //calculate products
    int inx_row = 0;
    for(int i = 0; i < num_shapes; i++, inx_row += num_shapes) {
      double* matrix_row = matrix[i];
      int shape0_inx = shape_inx[i];
      for(int k = 0; k < num_shapes; k++) {
        int shape1_inx = shape_inx[k];

        double value = 0.0;
        for(int j = 0; j < num_gip_points; j++) {
          double gip_x = gip_points[j][H2D_GIP2D_X], gip_y = gip_points[j][H2D_GIP2D_Y];
          double value0 = shapeset.get_value(H2D_FN_VALUE, shape0_inx, gip_x, gip_y, 0);
          double value1 = shapeset.get_value(H2D_FN_VALUE, shape1_inx, gip_x, gip_y, 0);
          double dx0 = shapeset.get_value(H2D_FN_DX, shape0_inx, gip_x, gip_y, 0);
          double dx1 = shapeset.get_value(H2D_FN_DX, shape1_inx, gip_x, gip_y, 0);
          double dy0 = shapeset.get_value(H2D_FN_DY, shape0_inx, gip_x, gip_y, 0);
          double dy1 = shapeset.get_value(H2D_FN_DY, shape1_inx, gip_x, gip_y, 0);

          value += gip_points[j][H2D_GIP2D_W] * (value0*value1 + dx0*dx1 + dy0*dy1);
        }

        matrix_row[k] = value;
      }
    }

    return matrix;
  }

  scalar H1ProjBasedSelector::evaluate_rsh_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx) {
    scalar total_value = 0;
    for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) {
      //get location and transform it
      double3 &gip_pt = sub_gip.gip_points[gip_inx];
      double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
      double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];

      //get value of a shape function
      scalar shape_value[3] = {0, 0, 0};
      shape_value[H2D_FN_VALUE] = shapeset->get_fn_value(shape_inx, ref_x, ref_y, 0);
      shape_value[H2D_FN_DX] = shapeset->get_dx_value(shape_inx, ref_x, ref_y, 0);
      shape_value[H2D_FN_DY] = shapeset->get_dy_value(shape_inx, ref_x, ref_y, 0);

      //get value of ref. solution
      scalar ref_value[3];
      ref_value[H2D_FN_VALUE] = sub_gip.rvals[H2D_FN_VALUE][gip_inx];
      ref_value[H2D_FN_DX] = sub_trf.coef_mx * sub_gip.rvals[H2D_FN_DX][gip_inx];
      ref_value[H2D_FN_DY] = sub_trf.coef_my * sub_gip.rvals[H2D_FN_DY][gip_inx];

      //evaluate a right-hand value
      scalar value = (shape_value[H2D_FN_VALUE] * ref_value[H2D_FN_VALUE])
        + (shape_value[H2D_FN_DX] * ref_value[H2D_FN_DX])
        + (shape_value[H2D_FN_DY] * ref_value[H2D_FN_DY]);

      total_value += gip_pt[H2D_GIP2D_W] * value;
    }
    return total_value;
  }

  double H1ProjBasedSelector::evaluate_error_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) {
    double total_error = 0;
    for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) {
      //get location and transform it
      double3 &gip_pt = sub_gip.gip_points[gip_inx];
      double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
      double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];

      //calculate value of projected solution
      scalar proj_value[3] = {0, 0, 0};
      for(int i = 0; i < elem_proj.num_shapes; i++) {
        int shape_inx = elem_proj.shape_inxs[i];
        proj_value[H2D_FN_VALUE] += elem_proj.shape_coefs[i] * shapeset->get_fn_value(shape_inx, ref_x, ref_y, 0);
        proj_value[H2D_FN_DX] += elem_proj.shape_coefs[i] * shapeset->get_dx_value(shape_inx, ref_x, ref_y, 0);
        proj_value[H2D_FN_DY] += elem_proj.shape_coefs[i] * shapeset->get_dy_value(shape_inx, ref_x, ref_y, 0);
      }

      //get value of ref. solution
      scalar ref_value[3];
      ref_value[H2D_FN_VALUE] = sub_gip.rvals[H2D_FN_VALUE][gip_inx];
      ref_value[H2D_FN_DX] = sub_trf.coef_mx * sub_gip.rvals[H2D_FN_DX][gip_inx];
      ref_value[H2D_FN_DY] = sub_trf.coef_my * sub_gip.rvals[H2D_FN_DY][gip_inx];

      //evaluate error
      double error = sqr(proj_value[H2D_FN_VALUE] - ref_value[H2D_FN_VALUE])
        + sqr(proj_value[H2D_FN_DX] - ref_value[H2D_FN_DX])
        + sqr(proj_value[H2D_FN_DY] - ref_value[H2D_FN_DY]);

      total_error += gip_pt[H2D_GIP2D_W] * error;
    }
    return total_error;
  }
}

