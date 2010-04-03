#include "../common.h"
#include "../matrix.h"
#include "../solution.h"
#include "../shapeset_hc_all.h"
#include "../element_to_refine.h"
#include "hcurl_proj_based_selector.h"

#define H2D_RV_VALUE0 0
#define H2D_RV_VALUE1 1
#define H2D_RV_D1DX 2
#define H2D_RV_D0DY 3

#ifdef COMPLEX

namespace RefinementSelectors {
  HcurlShapeset HcurlProjBasedSelector::default_shapeset;

  HcurlProjBasedSelector::HcurlProjBasedSelector(AdaptType adapt_type, double conv_exp, int max_order, HcurlShapeset* user_shapeset)
    : ProjBasedSelector(adapt_type, conv_exp, max_order, user_shapeset == NULL ? &default_shapeset : user_shapeset, 0) {}

  scalar** HcurlProjBasedSelector::precalc_ref_solution(int inx_son, Solution* rsln, Element* element, int intr_gip_order) {
    //set element and integration order
    rsln->set_active_element(element);
    rsln->set_quad_order(intr_gip_order);

    //fill with values
    scalar** rvals_son = precalc_rvals[inx_son];
    rvals_son[H2D_RV_VALUE0] = rsln->get_fn_values(0);
    rvals_son[H2D_RV_VALUE1] = rsln->get_fn_values(1);
    rvals_son[H2D_RV_D1DX] = rsln->get_dx_values(1);
    rvals_son[H2D_RV_D0DY] = rsln->get_dy_values(0);

    return rvals_son;
  }

  double** HcurlProjBasedSelector::build_projection_matrix(Shapeset& shapeset,
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
          double value0[2] = { shapeset.get_value(H2D_FN_VALUE, shape0_inx, gip_x, gip_y, 0), shapeset.get_value(H2D_FN_VALUE, shape0_inx, gip_x, gip_y, 1) };
          double value1[2] = { shapeset.get_value(H2D_FN_VALUE, shape1_inx, gip_x, gip_y, 0), shapeset.get_value(H2D_FN_VALUE, shape1_inx, gip_x, gip_y, 1) };
          double d1dx0 = shapeset.get_value(H2D_FN_DX, shape0_inx, gip_x, gip_y, 1);
          double d1dx1 = shapeset.get_value(H2D_FN_DX, shape1_inx, gip_x, gip_y, 1);
          double d0dy0 = shapeset.get_value(H2D_FN_DY, shape0_inx, gip_x, gip_y, 0);
          double d0dy1 = shapeset.get_value(H2D_FN_DY, shape1_inx, gip_x, gip_y, 0);
          double curl0 = d1dx0 - d0dy0;
          double curl1 = d1dx1 - d0dy1;

          value += gip_points[j][H2D_GIP2D_W] * (value0[0]*value1[0] + value0[1]*value1[1] + curl0*curl1);
        }

        matrix_row[k] = value;
      }
    }

    return matrix;
  }

  scalar HcurlProjBasedSelector::evaluate_rsh_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx) {
    double coef_curl = std::abs(sub_trf.coef_mx * sub_trf.coef_my);
    scalar total_value = 0;
    for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) {
      //get location and transform it
      double3 &gip_pt = sub_gip.gip_points[gip_inx];
      double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
      double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];

      //get value of a shape function
      scalar shape_value0 = shapeset->get_fn_value(shape_inx, ref_x, ref_y, 0);
      scalar shape_value1 = shapeset->get_fn_value(shape_inx, ref_x, ref_y, 1);
      scalar shape_curl = shapeset->get_dx_value(shape_inx, ref_x, ref_y, 1) - shapeset->get_dy_value(shape_inx, ref_x, ref_y, 0);

      //get value of ref. solution
      scalar ref_value0 = sub_trf.coef_mx * sub_gip.rvals[H2D_RV_VALUE0][gip_inx];
      scalar ref_value1 = sub_trf.coef_my * sub_gip.rvals[H2D_RV_VALUE1][gip_inx];
      scalar ref_curl = coef_curl * (sub_gip.rvals[H2D_RV_D1DX][gip_inx] - sub_gip.rvals[H2D_RV_D0DY][gip_inx]); //coef_curl * curl

      //evaluate a right-hand value
      scalar value = (shape_value0 * ref_value0)
        + (shape_value1 * ref_value1)
        + (shape_curl * ref_curl);

      total_value += gip_pt[H2D_GIP2D_W] * value;
    }
    return total_value;
  }

  double HcurlProjBasedSelector::evaluate_error_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) {
    double total_error = 0;
    double coef_curl = std::abs(sub_trf.coef_mx * sub_trf.coef_my);
    for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) {
      //get location and transform it
      double3 &gip_pt = sub_gip.gip_points[gip_inx];
      double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
      double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];

      //calculate value of projected solution
      scalar proj_value0 = 0, proj_value1 = 0, proj_curl = 0;
      for(int i = 0; i < elem_proj.num_shapes; i++) {
        int shape_inx = elem_proj.shape_inxs[i];
        proj_value0 += elem_proj.shape_coefs[i] * shapeset->get_fn_value(shape_inx, ref_x, ref_y, 0);
        proj_value1 += elem_proj.shape_coefs[i] * shapeset->get_fn_value(shape_inx, ref_x, ref_y, 1);
        proj_curl += elem_proj.shape_coefs[i] * (shapeset->get_dx_value(shape_inx, ref_x, ref_y, 1) - shapeset->get_dy_value(shape_inx, ref_x, ref_y, 0));
      }

      //get value of ref. solution
      scalar ref_value0 = sub_trf.coef_mx * sub_gip.rvals[H2D_RV_VALUE0][gip_inx];
      scalar ref_value1 = sub_trf.coef_my * sub_gip.rvals[H2D_RV_VALUE1][gip_inx];
      scalar ref_curl = coef_curl * (sub_gip.rvals[H2D_RV_D1DX][gip_inx] - sub_gip.rvals[H2D_RV_D0DY][gip_inx]); //coef_curl * curl

      //evaluate error
      double error = sqr(proj_value0 - ref_value0)
        + sqr(proj_value1 - ref_value1)
        + sqr(proj_curl - ref_curl);

      total_error += gip_pt[H2D_GIP2D_W] * error;
    }
    return total_error;
  }
}

#endif
