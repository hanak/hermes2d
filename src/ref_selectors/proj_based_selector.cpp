#include "../common.h"
#include "../solution.h"
#include "../linsystem.h"
#include "../quad_all.h"
#include "../element_to_refine.h"
#include "proj_based_selector.h"

namespace RefinementSelectors {
  
  ProjBasedSelector::ProjBasedSelector(AdaptType adapt_type, double conv_exp, int max_order, Shapeset* shapeset)
    : OptimumSelector(adapt_type, conv_exp, max_order, shapeset)
    , rhs_cache(NULL)
  {
    //build shape indices
    build_shape_indices(MODE_TRIANGLE);
    build_shape_indices(MODE_QUAD);

    //clear matrix cache
    for(int i = 0; i < H2DRS_MAX_ORDER+1; i++)
      for(int k = 0; k < H2DRS_MAX_ORDER+1; k++)
        proj_matrices[i][k] = NULL;

    //allocate cache
    int max_inx = max_shape_inx[0];
    for(int i = 1; i < H2D_NUM_MODES; i++)
      max_inx = std::max(max_inx, max_shape_inx[i]);
    rhs_cache = new ValueCacheItem<scalar>[max_inx + 1];
  }

  ProjBasedSelector::~ProjBasedSelector() {
    delete[] rhs_cache;
    for(int i = 0; i < H2DRS_MAX_ORDER+1; i++) {
      for(int k = 0; k < H2DRS_MAX_ORDER+1; k++) {
        if (proj_matrices[i][k] != NULL)
          free(proj_matrices[i][k]);

      }
    }
  }

  void ProjBasedSelector::build_shape_indices(const int mode) {
    std::vector<ShapeInx> &indices = shape_indices[mode];
    int* next_order = this->next_order_shape[mode];
    int& max_shape_inx = this->max_shape_inx[mode];
    int num_edges = (mode == MODE_QUAD) ? 4 : 3;
    shapeset->set_mode(mode);

    //cleanup
    indices.clear();
    indices.reserve((H2DRS_MAX_ORDER+1) * (H2DRS_MAX_ORDER+1));

    //TODO: what happens in a shapeset that lacks vertex-function (L2)
    //order 1
    next_order[0] = 0;
    for (int i = 0; i < num_edges; i++)
      indices.push_back(ShapeInx(1, 1, shapeset->get_vertex_index(i)));
    next_order[1] = (int)indices.size();

    //TODO: what happens in a shapeset that lacks edge-functions (HCurl)
    //order > 1
    max_shape_inx = 0;
    int examined_shape = 0;
    for (int i = 2; i <= H2DRS_MAX_ORDER; i++) {
      //edge functions
      if (mode == MODE_QUAD) {
        for (int j = 0; j < num_edges; j++)
          indices.push_back(ShapeInx(((j&1)==0) ? i : 0, ((j&1)!=0) ? i : 0, shapeset->get_edge_index(j, 0, i)));
      }
      else {
        for (int j = 0; j < num_edges; j++)
          indices.push_back(ShapeInx(i, i, shapeset->get_edge_index(j, 0, i)));
      }

      //bubble functions
      int bubble_order = (mode == MODE_QUAD) ? H2D_MAKE_QUAD_ORDER(i, i) : i;
      int num_bubbles = shapeset->get_num_bubbles(bubble_order);
      int* bubble_inxs = shapeset->get_bubble_indices(bubble_order);
      for(int j = 0; j < num_bubbles; j++) {
        int quad_order = shapeset->get_order(bubble_inxs[j]);
        int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
        if (std::max(order_h, order_v) == i)
          indices.push_back(ShapeInx(order_h, order_v, bubble_inxs[j]));
      }

      //store index of the next order
      next_order[i] = (int)indices.size();

      //update maximum
      while(examined_shape < next_order[i]) {
        max_shape_inx = std::max(max_shape_inx, indices[examined_shape].inx);
        examined_shape++;
      }
    }
  }

  void ProjBasedSelector::evaluate_cands_error(Element* e, Solution* rsln, double* avg_error, double* dev_error) {
    bool tri = e->is_triangle();

    // find range of orders
    int max_quad_order_h, max_quad_order_p, max_quad_order_aniso;
    evaluate_cands(&max_quad_order_h, &max_quad_order_p, &max_quad_order_aniso);

    // calculate (partial) projection errors
    SonProjectionError herr[4], anisoerr[4], perr;
    calc_projection_errors(e, max_quad_order_h, max_quad_order_p, max_quad_order_aniso, rsln, herr, anisoerr, perr);

    //evaluate errors and dofs
    double sum_err = 0.0;
    double sum_sqr_err = 0.0;
    int num_processed = 0;
    Cand& unrefined_c = candidates[0];
    for (unsigned i = 0; i < candidates.size(); i++) {
      Cand& c = candidates[i];
      if (tri) { //triangle
        switch(c.split) {
        case H2D_REFINEMENT_H:
          c.error = 0.0;
          for (int j = 0; j < H2D_MAX_ELEMENT_SONS; j++) {
            int order = H2D_GET_H_ORDER(c.p[j]);
            c.error += herr[j][order][order];
          }
          c.error *= 0.25; //convert to reference domain of a candidate because every son of a candidate is integrated over a reference domain;
          break;

        case H2D_REFINEMENT_P:
          {
            int order = H2D_GET_H_ORDER(c.p[0]);
            c.error = perr[order][order];
          }
          break;

        default:
          error("Unknown split type \"%d\" at candidate %d", c.split, i);
        }
      }
      else { //quad
        switch(c.split) {
        case H2D_REFINEMENT_H:
          c.error = 0.0;
          for (int j = 0; j < H2D_MAX_ELEMENT_SONS; j++) {
            int order_h = H2D_GET_H_ORDER(c.p[j]), order_v = H2D_GET_V_ORDER(c.p[j]);
            c.error += herr[j][order_h][order_v];
          }
          c.error *= 0.25; //convert to reference domain of a candidate because every son of a candidate is integrated over a reference domain;
          break;

        case H2D_REFINEMENT_ANISO_H:
        case H2D_REFINEMENT_ANISO_V:
          {
            c.error = 0.0;
            for (int j = 0; j < 2; j++)
              c.error += anisoerr[(c.split == H2D_REFINEMENT_ANISO_H) ? j : j+2][H2D_GET_H_ORDER(c.p[j])][H2D_GET_V_ORDER(c.p[j])];
            c.error *= 0.5;  //convert to reference domain of a candidate because every son of a candidate is integrated over a reference domain;
          }
          break;

        case H2D_REFINEMENT_P:
          {
            int order_h = H2D_GET_H_ORDER(c.p[0]), order_v = H2D_GET_V_ORDER(c.p[0]);
            c.error = perr[order_h][order_v];
          }
          break;

        default:
          error("Unknown split type \"%d\" at candidate %d", c.split, i);
        }
      }

      c.error = sqrt(c.error);
      if (!i || c.error <= unrefined_c.error)
      {
        sum_err += log10(c.error);
        sum_sqr_err += sqr(log10(c.error));
        num_processed++;
      }
    }

    *avg_error = sum_err / num_processed;  // mean
    *dev_error = sqrt(sum_sqr_err/num_processed - sqr(*avg_error)); // deviation is square root of variance
  }

  void ProjBasedSelector::calc_projection_errors(Element* e, const int max_quad_order_h, const int max_quad_order_p, const int max_quad_order_aniso, Solution* rsln, SonProjectionError herr[4], SonProjectionError anisoerr[4], SonProjectionError perr) {
    assert_msg(H2D_GET_H_ORDER(max_quad_order_h) <= H2DRS_MAX_ORDER && H2D_GET_V_ORDER(max_quad_order_h) <= H2DRS_MAX_ORDER, "E maximum allowed order of a son of H-candidate is %d but order (H:%d,V:%d) requested", H2DRS_MAX_ORDER, H2D_GET_H_ORDER(max_quad_order_h), H2D_GET_V_ORDER(max_quad_order_h));
    assert_msg(H2D_GET_H_ORDER(max_quad_order_p) <= H2DRS_MAX_ORDER && H2D_GET_V_ORDER(max_quad_order_p) <= H2DRS_MAX_ORDER, "E maximum allowed order of a son of P-candidate is %d but order (H:%d,V:%d) requested", H2DRS_MAX_ORDER, H2D_GET_H_ORDER(max_quad_order_p), H2D_GET_V_ORDER(max_quad_order_p));
    assert_msg(H2D_GET_H_ORDER(max_quad_order_aniso) <= H2DRS_MAX_ORDER && H2D_GET_V_ORDER(max_quad_order_aniso) <= H2DRS_MAX_ORDER, "E maximum allowed order of a son of ANISO-candidate is %d but order (H:%d,V:%d) requested", H2DRS_MAX_ORDER, H2D_GET_H_ORDER(max_quad_order_aniso), H2D_GET_V_ORDER(max_quad_order_aniso));

    int mode = e->get_mode();

    // select quadrature, obtain integration points and weights
    Quad2D* quad = &g_quad_2d_std;
    quad->set_mode(mode);
    rsln->set_quad_2d(quad);
    //TODO: select GIP order according to really needed order
    double3* gip_points = quad->get_points(H2DRS_INTR_GIP_ORDER);
    int num_gip_points = quad->get_num_points(H2DRS_INTR_GIP_ORDER);

    // everything is done on the reference domain
    rsln->enable_transform(false);

    // obtain reference solution values on all four refined sons
    scalar* rval[H2D_MAX_ELEMENT_SONS][3];
    Element* base_element = rsln->get_mesh()->get_element(e->id);
    assert(!base_element->active);
    for (int son = 0; son < H2D_MAX_ELEMENT_SONS; son++)
    {
      Element* e = base_element->sons[son];
      assert(e != NULL);
      rsln->set_active_element(e);
      rsln->set_quad_order(H2DRS_INTR_GIP_ORDER);
      rval[son][H2D_FN_VALUE] = rsln->get_fn_values();
      rval[son][H2D_FN_DX] = rsln->get_dx_values();
      rval[son][H2D_FN_DY] = rsln->get_dy_values();
    }

    //H-candidates
    if (max_quad_order_h != 0) {
      Trf trf_identity = { {1.0, 1.0}, {0.0, 0.0} };
      Trf* p_trf_identity[1] = { &trf_identity };
      double coef_mm = 1;
      for(int son = 0; son < H2D_MAX_ELEMENT_SONS; son++) {
        scalar **sub_rval[1] = { rval[son] };
        calc_proj_error_cand_son(mode, gip_points, num_gip_points
          , 1, &base_element->sons[son], p_trf_identity, sub_rval, &coef_mm, &coef_mm
          , max_quad_order_h, herr[son]);
      }
    }

    //ANISO-candidates
    if (max_quad_order_aniso != 0) {
      const double mx[4] = { 2.0, 2.0, 1.0, 1.0}; //scale coefficients of dx for X-axis due to trasformations
      const double my[4] = { 1.0, 1.0, 2.0, 2.0}; //scale coefficients of dy for Y-axis due to trasformations
      const int sons[4][2] = { {0,1}, {3,2}, {0,3}, {1,2} }; //indices of sons for sub-areas
      const int tr[4][2]   = { {6,7}, {6,7}, {4,5}, {4,5} }; //indices of ref. domain transformations for sub-areas
      for(int version = 0; version < 4; version++) { // 2 sons for vertical split, 2 sons for horizontal split
        Trf* sub_trfs[2] = { &quad_trf[tr[version][0]], &quad_trf[tr[version][1]] };
        Element* sub_elems[2] = { base_element->sons[sons[version][0]], base_element->sons[sons[version][1]] };
        scalar **sub_rval[2] = { rval[sons[version][0]], rval[sons[version][1]] };
        double coefs_mx[2] = { mx[version], mx[version] }, coefs_my[2] = { my[version], my[version] };
        calc_proj_error_cand_son(mode, gip_points, num_gip_points
          , 2, sub_elems, sub_trfs, sub_rval, coefs_mx, coefs_my
          , max_quad_order_aniso, anisoerr[version]);
      }
    }

    //P-candidates
    if (max_quad_order_p != 0) {
      Trf* src_trfs = NULL;
      if (mode == MODE_TRIANGLE)
        src_trfs = tri_trf;
      else
        src_trfs = quad_trf;
      Trf* sub_trfs[4] = { &src_trfs[0], &src_trfs[1], &src_trfs[2], &src_trfs[3] };
      scalar **sub_rval[4] = { rval[0], rval[1], rval[2], rval[3] };
      double coefs_mm[4] = { 2.0, 2.0, 2.0, (mode == MODE_TRIANGLE) ? -2.0 : 2.0 };
      calc_proj_error_cand_son(mode, gip_points, num_gip_points
        , 4, base_element->sons, sub_trfs, sub_rval, coefs_mm, coefs_mm
        , max_quad_order_p, perr);
    }
  }

  void ProjBasedSelector::calc_proj_error_cand_son(const int mode
    , double3* gip_points, int num_gip_points
    , const int num_sub, Element** sub_elems, Trf** sub_trfs, scalar*** sub_rvals, double* coefs_mx, double* coefs_my
    , int max_quad_order
    , SonProjectionError errors
    ) {
    //allocate space
    int max_num_shapes = next_order_shape[mode][current_max_order] - next_order_shape[mode][0];
    scalar* right_side = new scalar[max_num_shapes];
    int* shape_inxs = new int[max_num_shapes];
    int* indx = new int[max_num_shapes]; //solver data
    double* d = new double[max_num_shapes]; //solver data
    double** proj_matrix = new_matrix<double>(max_num_shapes, max_num_shapes);
    std::vector<ShapeInx>& full_shape_indices = shape_indices[mode];

    //clenup of the cache
    memset(rhs_cache, 0, sizeof(ValueCacheItem<scalar>) * max_shape_inx[mode]); //TODO: can be done efficient by remembering the previous order during creating a list of indices for current range of orders

    //calculate for all orders
    double sub_area_corr_coef = 1.0 / num_sub;
    OrderPermutator order_perm(H2D_MAKE_QUAD_ORDER(1, 1), max_quad_order, mode == MODE_TRIANGLE);
    do {
      int quad_order = order_perm.get_quad_order();
      int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);

      //build a list of shape indices from the full list
      int num_shapes = 0;
      unsigned int inx_shape = 0;
      while (inx_shape < full_shape_indices.size()) {
        ShapeInx& shape = full_shape_indices[inx_shape];
        if (order_h >= shape.order_h && order_v >= shape.order_v) {
          debug_assert(num_shapes < max_num_shapes, "E more shapes than predicted, possible incosistency");
          shape_inxs[num_shapes] = shape.inx;
          num_shapes++;
        }
        inx_shape++;
      }

      //calculate projection matrix
      if (proj_matrices[order_h][order_v] == NULL)
        proj_matrices[order_h][order_v] = build_projection_matrix(*shapeset, gip_points, num_gip_points, shape_inxs, num_shapes);
      copy_matrix(proj_matrix, proj_matrices[order_h][order_v], num_shapes, num_shapes); //copy projection matrix because original matrix will be modified

      //build right side (fill cache values that are missing)
      for(int inx_sub = 0; inx_sub < num_sub; inx_sub++) {
        Element* sub_elem = sub_elems[inx_sub];
        ElemSubTrf sub_trf = { sub_trfs[inx_sub], coefs_mx[inx_sub], coefs_my[inx_sub] };
        ElemGIP sub_gip = { gip_points, num_gip_points, sub_rvals[inx_sub] };

        for(int k = 0; k < num_shapes; k++) {
          int shape_inx = shape_inxs[k];
          ValueCacheItem<scalar>& shape_rhs_cache = rhs_cache[shape_inx];
          if (!shape_rhs_cache.is_valid())
            shape_rhs_cache.set(shape_rhs_cache.get() + evaluate_rsh_sub_element(sub_elem, sub_gip, sub_trf, shape_inx));
        }
      }

      //copy values from cache and apply area correction coefficient
      for(int k = 0; k < num_shapes; k++) {
        ValueCacheItem<scalar>& shape_rhs_cache = rhs_cache[shape_inxs[k]];
        right_side[k] = sub_area_corr_coef * shape_rhs_cache.get();
        shape_rhs_cache.mark();
      }

      //solve
      ludcmp(proj_matrix, num_shapes, indx, d);
      lubksb<scalar>(proj_matrix, num_shapes, indx, right_side);

      //calculate error
      double error = 0;
      for(int inx_sub = 0; inx_sub < num_sub; inx_sub++) {
        Element* sub_elem = sub_elems[inx_sub];
        Trf* ref_coord_transf = sub_trfs[inx_sub];
        double coef_mx = coefs_mx[inx_sub], coef_my = coefs_my[inx_sub];
        ElemSubTrf sub_trf = { sub_trfs[inx_sub], coefs_mx[inx_sub], coefs_my[inx_sub] };
        ElemGIP sub_gip = { gip_points, num_gip_points, sub_rvals[inx_sub] };
        ElemProj elem_proj = { shape_inxs, num_shapes, right_side, quad_order };

        error += evaluate_error_sub_element(sub_elem, sub_gip, sub_trf, elem_proj);
      }
      errors[order_h][order_v] = error * sub_area_corr_coef; //apply area correction coefficient
    } while (order_perm.next());

    //clenaup
    free(proj_matrix);
    delete[] right_side;
    delete[] shape_inxs;
    delete[] indx;
    delete[] d;
  }

}
