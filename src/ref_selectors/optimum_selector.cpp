#include "../common.h"
#include "../solution.h"
#include "../linsystem.h"
#include "../quad_all.h"
#include "../element_to_refine.h"
#include "optimum_selector.h"

#define H2DST_ANY H2DST_VERTEX | H2DST_HORIZ_EDGE | H2DST_VERT_EDGE | H2DST_TRI_EDGE | H2DST_BUBBLE ///< Any type of shape. Used just for masky

namespace RefinementSelectors {

  HERMES2D_API std::ostream &operator<<(std::ostream &stream, const OptimumSelector::Cand& cand) {
    stream.precision(2);
    stream << "split:" << get_refin_str(cand.split);
    stream << " err:" << std::scientific << cand.error << " dofs:" << cand.dofs << " ";

    int num_sons = cand.get_num_sons();
    stream << "[";
    for(int i = 0; i < num_sons; i++) {
      if (i > 0)
        stream << " ";
      stream << get_quad_order_str(cand.p[i]);
    }
    stream << "]";
    return stream;
  }

  HERMES2D_API bool is_hp(const CandList cand_list) {
    switch(cand_list) {
      case H2D_P_ISO:
      case H2D_P_ANISO:
      case H2D_H_ISO:
      case H2D_H_ANISO: return false; break;
      case H2D_HP_ISO:
      case H2D_HP_ANISO_H:
      case H2D_HP_ANISO_P:
      case H2D_HP_ANISO: return true; break;
      default: error("invalid adapt type %d", cand_list); return false;
    }
  }

  OrderPermutator::OrderPermutator(int start_quad_order, int end_quad_order, bool iso_p, int* tgt_quad_order)
    : start_order_h(H2D_GET_H_ORDER(start_quad_order)), start_order_v(H2D_GET_V_ORDER(start_quad_order))
    , end_order_h(H2D_GET_H_ORDER(end_quad_order)), end_order_v(H2D_GET_V_ORDER(end_quad_order))
    , iso_p(iso_p), tgt_quad_order(tgt_quad_order) {
    assert_msg(start_order_h <= end_order_h && start_order_v <= end_order_v, "End orders (H:%d, V:%d) are below start orders (H:%d, V:%d).", end_order_h, end_order_v, start_order_h, start_order_v);
    reset();
  }

  bool OrderPermutator::next() {
    if (order_h >= end_order_h && order_v >= end_order_v)
      return false;
    else {
      if (iso_p) {
        order_h++; order_v++;
      }
      else {
        order_h++;
        if (order_h > end_order_h) {
          order_h = start_order_h;
          order_v++;
        }
      }

      if (tgt_quad_order != NULL)
        *tgt_quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
      return true;
    }
  }

  void OrderPermutator::reset() {
    order_h = start_order_h;
    order_v = start_order_v;
    if (tgt_quad_order != NULL)
      *tgt_quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
  }

  OptimumSelector::OptimumSelector(CandList cand_list, double conv_exp, int max_order, Shapeset* shapeset, const Range<int>& vertex_order, const Range<int>& edge_bubble_order)
    : Selector(max_order), cand_list(cand_list)
    , conv_exp(conv_exp), shapeset(shapeset) {
    error_if(shapeset == NULL, "Shapeset is NULL.");

    //build shape indices
    build_shape_indices(MODE_TRIANGLE, vertex_order, edge_bubble_order);
    build_shape_indices(MODE_QUAD, vertex_order, edge_bubble_order);
  }

  void OptimumSelector::build_shape_indices(const int mode, const Range<int>& vertex_order, const Range<int>& edge_bubble_order) {
    std::vector<ShapeInx> &indices = shape_indices[mode];
    int* next_order = this->next_order_shape[mode];
    int& max_shape_inx = this->max_shape_inx[mode];
    int num_edges = (mode == MODE_QUAD) ? 4 : 3;
    shapeset->set_mode(mode);
    bool &has_vertex = has_vertex_shape[mode];
    bool &has_edge = has_edge_shape[mode];
    bool &has_bubble = has_bubble_shape[mode];

    //cleanup
    indices.clear();
    indices.reserve((H2DRS_MAX_ORDER+1) * (H2DRS_MAX_ORDER+1));
    has_vertex = has_edge = has_bubble = false;

    //get total range of orders
    Range<int> order_range = Range<int>::make_envelope(vertex_order, edge_bubble_order);

    //for all orders
    max_shape_inx = 0;
    int examined_shape = 0;
    for(int i = order_range.lower(); i <= order_range.upper(); i++) {
      //vertex functions
      if (vertex_order.is_in_closed(i)) {
        for (int i = 0; i < num_edges; i++) {
          int inx = shapeset->get_vertex_index(i);
          if (inx >= 0) {
            indices.push_back(ShapeInx(1, 1, inx, H2DST_VERTEX));
            has_vertex = true;
          }
        }
      }

      //edge functions and bubble functions
      if (edge_bubble_order.is_in_closed(i)) {
        //edge functions
        if (mode == MODE_QUAD) {
          for (int j = 0; j < num_edges; j++) {
            int inx = shapeset->get_edge_index(j, 0, i);
            if (inx >= 0) {
              if ((j&1) == 0) //horizontal edge
                indices.push_back(ShapeInx(i, 0, inx, H2DST_HORIZ_EDGE));
              else //vertical edge
                indices.push_back(ShapeInx(0, i, inx, H2DST_VERT_EDGE));
              has_edge = true;
            }
          }
        }
        else {
          for (int j = 0; j < num_edges; j++) {
            int inx = shapeset->get_edge_index(j, 0, i);
            if (inx >= 0) {
              indices.push_back(ShapeInx(i, i, inx, H2DST_TRI_EDGE));
              has_edge = true;
            }
          }
        }

        //bubble functions
        int bubble_order = (mode == MODE_QUAD) ? H2D_MAKE_QUAD_ORDER(i, i) : i;
        int num_bubbles = shapeset->get_num_bubbles(bubble_order);
        int* bubble_inxs = shapeset->get_bubble_indices(bubble_order);
        for(int j = 0; j < num_bubbles; j++) {
          //retrieve order of a bubble function
          int quad_order = shapeset->get_order(bubble_inxs[j]);
          int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
          warn_if(std::max(order_h, order_v) > i, "requested bubble functions of up to order %d but retrieved bubble function %d has order (H:%d; V:%d), mode %d", i, bubble_inxs[j], order_h, order_v, mode);

          //HACK: check whether the bubble function is not already used (Hcurl does that)
          int inx_bubble = bubble_inxs[j];
          std::vector<ShapeInx>::const_iterator inx = indices.begin();
          while (inx != indices.end() && inx->inx != inx_bubble)
            inx++;
          if (inx == indices.end()) { //bubble is unique: add it
            //HACK: Hcurl returns bubble function 2/0 when bubble functions of up to order 1/1 are requested
            order_h = std::min(order_h, i);
            order_v = std::min(order_v, i);

            //add bubble function
            indices.push_back(ShapeInx(order_h, order_v, bubble_inxs[j], H2DST_BUBBLE));
            has_bubble = true;
          }
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
  }

  int OptimumSelector::calc_num_shapes(int mode, int order_h, int order_v, int allowed_type_mask) {
    //test whether the evaluation is necessary
    bool full_eval = false;
    if ((allowed_type_mask & H2DST_VERTEX) != 0)
      full_eval |= has_vertex_shape[mode];
    if ((allowed_type_mask & (H2DST_HORIZ_EDGE | H2DST_VERT_EDGE | H2DST_TRI_EDGE)) != 0)
      full_eval |= has_edge_shape[mode];
    if ((allowed_type_mask & H2DST_BUBBLE) != 0)
      full_eval |= has_bubble_shape[mode];

    //evaluate
    if (full_eval) {
      std::vector<ShapeInx>& shapes = shape_indices[mode];
      int num = 0;
      std::vector<ShapeInx>::const_iterator shape = shapes.begin();
      while (shape != shapes.end()) {
        if (((int)shape->type & allowed_type_mask) != 0) {
          if ((order_h == H2DRS_ORDER_ANY || shape->order_h <= order_h) && (order_v == H2DRS_ORDER_ANY || shape->order_v <= order_v))
            num++;
        }
        shape++;
      }
      return num;
    }
    else
      return 0;
  }

  void OptimumSelector::append_candidates_split(const int start_quad_order, const int last_quad_order, const int split, bool iso_p) {
    //check whether end orders are not lower than start orders
    if (last_quad_order < 0 || start_quad_order < 0)
      return;
    if (H2D_GET_H_ORDER(start_quad_order) > H2D_GET_H_ORDER(last_quad_order) || H2D_GET_V_ORDER(start_quad_order) > H2D_GET_V_ORDER(last_quad_order))
      return;

    //get number of sons
    const int num_sons = get_refin_sons(split);

    //initialize orders
    int quad_orders[H2D_MAX_ELEMENT_SONS];
    OrderPermutator quad_perms[H2D_MAX_ELEMENT_SONS];
    for(int i = 0; i < num_sons; i++) {
      quad_orders[i] = start_quad_order;
      quad_perms[i] = OrderPermutator(start_quad_order, last_quad_order, iso_p, &quad_orders[i]);
    }
    for(int i = num_sons; i < H2D_MAX_ELEMENT_SONS; i++)
      quad_orders[i] = 0;

    //generate permutations of orders
    bool quit = false;
    while(!quit) {
      do { //create permutation of son 0
        candidates.push_back(Cand(split, quad_orders));
      } while (quad_perms[0].next());

      //reset son 0
      quad_perms[0].reset();

      //increment orders of other sons
      int inx_son = 1;
      while (inx_son < num_sons && !quad_perms[inx_son].next()) {
        quad_perms[inx_son].reset(); //reset order of the son
        inx_son++;
      }
      if (inx_son >= num_sons)
        quit = true;
    }
  }

  void OptimumSelector::create_candidates(Element* e, int quad_order, int max_ha_quad_order, int max_p_quad_order) {
    int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
    int max_p_order_h = H2D_GET_H_ORDER(max_p_quad_order), max_p_order_v = H2D_GET_V_ORDER(max_p_quad_order);
    int max_ha_order_h = H2D_GET_H_ORDER(max_ha_quad_order), max_ha_order_v = H2D_GET_V_ORDER(max_ha_quad_order);
    bool tri = e->is_triangle();

    //clear list of candidates
    candidates.clear();
    if (candidates.capacity() < H2DRS_ASSUMED_MAX_CANDS)
      candidates.reserve(H2DRS_ASSUMED_MAX_CANDS);

    //generate all P-candidates (start from intention of generating all possible candidates and restrict it according to the given adapt-type)
    bool iso_p = false;
    int start_quad_order = quad_order;
    int last_quad_order = H2D_MAKE_QUAD_ORDER(std::min(max_p_order_h, order_h+H2DRS_MAX_ORDER_INC), std::min(max_p_order_v, order_v+H2DRS_MAX_ORDER_INC));
    switch(cand_list) {
      case H2D_H_ISO:
      case H2D_H_ANISO:
        last_quad_order = start_quad_order; break; //no P-candidates except the original candidate

      case H2D_P_ISO:
      case H2D_HP_ISO:
      case H2D_HP_ANISO_H:
        iso_p = true; break; //iso change of orders
    }
    append_candidates_split(quad_order, last_quad_order, H2D_REFINEMENT_P, tri || iso_p);

    //generate all H-candidates
    iso_p = false;
    int start_order_h = std::max(current_min_order, (order_h+1) / 2), start_order_v = std::max(current_min_order, (order_v+1) / 2);
    start_quad_order = H2D_MAKE_QUAD_ORDER(start_order_h, start_order_v);
    last_quad_order = H2D_MAKE_QUAD_ORDER(std::min(max_ha_order_h, std::min(start_order_h + H2DRS_MAX_ORDER_INC, order_h)), std::min(max_ha_order_v, std::min(start_order_v + H2DRS_MAX_ORDER_INC, order_v)));
    switch(cand_list) {
      case H2D_H_ISO:
      case H2D_H_ANISO:
        last_quad_order = start_quad_order; break; //no only one candidate will be created

      case H2D_P_ISO:
      case H2D_P_ANISO:
        last_quad_order = -1; break; //no H-candidate will be generated

      case H2D_HP_ISO:
      case H2D_HP_ANISO_H:
        iso_p = true; break; //iso change of orders
    }
    append_candidates_split(start_quad_order, last_quad_order, H2D_REFINEMENT_H, tri || iso_p);

    //generate all ANISO-candidates
    if (!tri && e->iro_cache < 8 //TODO: find out what iro_cache does and why 8
      && (cand_list == H2D_H_ANISO || cand_list == H2D_HP_ANISO_H || cand_list == H2D_HP_ANISO)) {
      iso_p = false;
      int start_quad_order_hz = H2D_MAKE_QUAD_ORDER(order_h, std::max(current_min_order, (order_v+1) / 2));
      int last_quad_order_hz = H2D_MAKE_QUAD_ORDER(std::min(max_ha_order_h, order_h+H2DRS_MAX_ORDER_INC), std::min(order_v, H2D_GET_V_ORDER(start_quad_order)+H2DRS_MAX_ORDER_INC));
      int start_quad_order_vt = H2D_MAKE_QUAD_ORDER(std::max(current_min_order, (order_h+1) / 2), order_v);
      int last_quad_order_vt = H2D_MAKE_QUAD_ORDER(std::min(order_h, H2D_GET_H_ORDER(start_quad_order)+H2DRS_MAX_ORDER_INC), std::min(max_ha_order_v, order_v+H2DRS_MAX_ORDER_INC));
      switch(cand_list) {
        case H2D_H_ANISO:
          last_quad_order_hz = start_quad_order_hz;
          last_quad_order_vt = start_quad_order_vt;
          break; //only one candidate will be created

        case H2D_HP_ANISO_H:
          iso_p = true; break; //iso change of orders
      }
      if (iso_p) { //make orders uniform: take mininmum order since nonuniformity is caused by different handling of orders along directions
        int order = std::min(H2D_GET_H_ORDER(start_quad_order_hz), H2D_GET_V_ORDER(start_quad_order_hz));
        start_quad_order_hz = H2D_MAKE_QUAD_ORDER(order, order);
        order = std::min(H2D_GET_H_ORDER(start_quad_order_vt), H2D_GET_V_ORDER(start_quad_order_vt));
        start_quad_order_vt = H2D_MAKE_QUAD_ORDER(order, order);

        order = std::min(H2D_GET_H_ORDER(last_quad_order_hz), H2D_GET_V_ORDER(last_quad_order_hz));
        last_quad_order_hz = H2D_MAKE_QUAD_ORDER(order, order);
        order = std::min(H2D_GET_H_ORDER(last_quad_order_vt), H2D_GET_V_ORDER(last_quad_order_vt));
        last_quad_order_vt = H2D_MAKE_QUAD_ORDER(order, order);
      }
      append_candidates_split(start_quad_order_hz, last_quad_order_hz, H2D_REFINEMENT_ANISO_H, iso_p);
      append_candidates_split(start_quad_order_vt, last_quad_order_vt, H2D_REFINEMENT_ANISO_V, iso_p);
    }
  }

  void OptimumSelector::update_cands_info(CandsInfo& info_h, CandsInfo& info_p, CandsInfo& info_aniso) const {
    std::vector<Cand>::const_iterator cand = candidates.begin();
    while (cand != candidates.end()) {
      CandsInfo* info = NULL;
      if (cand->split == H2D_REFINEMENT_H) info = &info_h;
      else if (cand->split == H2D_REFINEMENT_P) info = &info_p;
      else if (cand->split == H2D_REFINEMENT_ANISO_H || cand->split == H2D_REFINEMENT_ANISO_V) info = &info_aniso;
      else { error("Invalid candidate type: %d.", cand->split); };

      //evaluate sons of candidates
      const int num_sons = cand->get_num_sons();
      for(int i = 0; i < num_sons; i++) {
        int son_order_h = H2D_GET_H_ORDER(cand->p[i]), son_order_v = H2D_GET_H_ORDER(cand->p[i]);
        if (son_order_h != son_order_v)
          info->uniform_orders = false;
        if (info->min_quad_order < 0 || info->max_quad_order < 0)
          info->min_quad_order = info->max_quad_order = H2D_MAKE_QUAD_ORDER(son_order_h, son_order_v);
        else {
          info->min_quad_order = H2D_MAKE_QUAD_ORDER(std::min(H2D_GET_H_ORDER(info->min_quad_order), son_order_h), std::min(H2D_GET_V_ORDER(info->min_quad_order), son_order_v));
          info->max_quad_order = H2D_MAKE_QUAD_ORDER(std::max(H2D_GET_H_ORDER(info->max_quad_order), son_order_h), std::max(H2D_GET_V_ORDER(info->max_quad_order), son_order_v));
        }
      }

      //next candidate
      cand++;
    }
  }

  void OptimumSelector::evaluate_cands_dof(Element* e, Solution* rsln) {
    bool tri = e->is_triangle();

    for (unsigned i = 0; i < candidates.size(); i++) {
      Cand& c = candidates[i];
      if (tri) { //triangle
        switch(c.split) {
        case H2D_REFINEMENT_H:
          {
            int central = 1;
            c.dofs = 0;
            for(int j = 0; j < H2D_MAX_ELEMENT_SONS; j++) {
              c.dofs += calc_num_shapes(MODE_TRIANGLE, H2D_GET_H_ORDER(c.p[j]), H2DRS_ORDER_ANY, H2DST_ANY);
              if (j != central)
                c.dofs -= calc_num_shapes(MODE_TRIANGLE, std::min(H2D_GET_H_ORDER(c.p[j]), H2D_GET_H_ORDER(c.p[central])), H2DRS_ORDER_ANY, H2DST_TRI_EDGE) / 3; //shared edge: since triangle has three edges which are identified by a single order this will find 3 x different edge of a given order
            }
            if (has_vertex_shape[MODE_TRIANGLE])
              c.dofs -= 2*3; //every vertex functions is added 3-times
          }
          break;

        case H2D_REFINEMENT_P:
          c.dofs = calc_num_shapes(MODE_TRIANGLE, H2D_GET_H_ORDER(c.p[0]), H2DRS_ORDER_ANY, H2DST_ANY);
          break;

        default:
          error("Unknown split type \"%d\" at candidate %d", c.split, i);
        }
      }
      else { //quad
        switch(c.split) {
        case H2D_REFINEMENT_H:
          c.dofs = 0;
          for(int j = 0; j < H2D_MAX_ELEMENT_SONS; j++)
            c.dofs += calc_num_shapes(MODE_QUAD, H2D_GET_H_ORDER(c.p[j]), H2D_GET_V_ORDER(c.p[j]), H2DST_ANY);
          for(int j = 0; j < 2; j++) { //shared edge functions
            c.dofs -= calc_num_shapes(MODE_QUAD, H2DRS_ORDER_ANY, std::min(H2D_GET_V_ORDER(c.p[2*j]), H2D_GET_V_ORDER(c.p[2*j + 1])), H2DST_VERT_EDGE) / 2; //shared vertical edge functions: every edge is twice there
            c.dofs -= calc_num_shapes(MODE_QUAD, std::min(H2D_GET_H_ORDER(c.p[j]), H2D_GET_H_ORDER(c.p[j^3])), H2DRS_ORDER_ANY, H2DST_HORIZ_EDGE) / 2; //shared horizontal edge functions: every edge is twice there
          }
          if (has_vertex_shape[MODE_QUAD])
            c.dofs -= 4 + 3; //edge vertex + central vertex
          break;

        case H2D_REFINEMENT_ANISO_H:
          c.dofs = calc_num_shapes(MODE_QUAD, H2D_GET_H_ORDER(c.p[0]), H2D_GET_V_ORDER(c.p[0]), H2DST_ANY);
          c.dofs += calc_num_shapes(MODE_QUAD, H2D_GET_H_ORDER(c.p[1]), H2D_GET_V_ORDER(c.p[1]), H2DST_ANY);
          c.dofs -= calc_num_shapes(MODE_QUAD, std::min(H2D_GET_H_ORDER(c.p[0]), H2D_GET_H_ORDER(c.p[1])), H2DRS_ORDER_ANY, H2DST_HORIZ_EDGE) / 2; //shared edge functions
          if (has_vertex_shape[MODE_QUAD])
            c.dofs -= 2; //shared vertex functions
          break;

        case H2D_REFINEMENT_ANISO_V:
          c.dofs = calc_num_shapes(MODE_QUAD, H2D_GET_H_ORDER(c.p[0]), H2D_GET_V_ORDER(c.p[0]), H2DST_ANY);
          c.dofs += calc_num_shapes(MODE_QUAD, H2D_GET_H_ORDER(c.p[1]), H2D_GET_V_ORDER(c.p[1]), H2DST_ANY);
          c.dofs -= calc_num_shapes(MODE_QUAD, H2DRS_ORDER_ANY, std::min(H2D_GET_V_ORDER(c.p[0]), H2D_GET_V_ORDER(c.p[1])), H2DST_VERT_EDGE) / 2; //shared edge functions
          if (has_vertex_shape[MODE_QUAD])
            c.dofs -= 2; //shared vertex functions
          break;

        case H2D_REFINEMENT_P:
          c.dofs = calc_num_shapes(MODE_QUAD, H2D_GET_H_ORDER(c.p[0]), H2D_GET_V_ORDER(c.p[0]), H2DST_ANY);
          break;

        default:
          error("Unknown split type \"%d\" at candidate %d", c.split, i);
        }
      }
    }
  }

  void OptimumSelector::evaluate_candidates(Element* e, Solution* rsln, double* avg_error, double* dev_error) {
    evaluate_cands_error(e, rsln, avg_error, dev_error);
    evaluate_cands_dof(e, rsln);
  }

  bool OptimumSelector::compare_cand_score(const Cand& a, const Cand& b) {
    return a.score > b.score;
  }

  void OptimumSelector::select_best_candidate(Element* e, const double avg_error, const double dev_error, int* selected_cand, int* selected_h_cand) {
    // select an above-average candidate with the steepest error decrease

    //calculate score of candidates
    Cand& unrefined = candidates[0];
    const int num_cands = (int)candidates.size();
    unrefined.score = 0;
    for (int i = 1; i < num_cands; i++) {
      Cand& cand = candidates[i];
      if (cand.error < unrefined.error && cand.dofs > unrefined.dofs)
        candidates[i].score = (log10(unrefined.error) - log10(cand.error)) / pow(cand.dofs - unrefined.dofs, conv_exp);
        //candidates[i].score = (log10(unrefined.error) - log10(cand.error)) / (pow(cand.dofs, conv_exp) - pow(unrefined.dofs, conv_exp)); //DEBUG: fixed
      else
        candidates[i].score = 0;
    }

    //sort according to the score
    if (num_cands > 2)
      std::sort(candidates.begin()+1, candidates.end(), compare_cand_score);

    //find first valid score that diffres from the next scores
    int imax = 1;
    while ((imax+1) < num_cands && std::abs(candidates[imax].score - candidates[imax+1].score) < H2DRS_SCORE_DIFF_ZERO) {
      //find the first candidate with a different score
      Cand& cand_current = candidates[imax];
      int imax_end = imax + 2;
      while (imax_end < num_cands && std::abs(cand_current.score - candidates[imax_end].score) < H2DRS_SCORE_DIFF_ZERO)
        imax_end++;

      imax = imax_end;
    }

    //find valid H-refinement candidate
    int h_imax = imax;
    while (h_imax < num_cands && candidates[h_imax].split != H2D_REFINEMENT_H)
      h_imax++;

    //make sure the result is valid: index is valid, selected candidate has a valid score
    if (imax >= num_cands || candidates[imax].score == 0)
      imax = 0;
    if (h_imax >= num_cands || candidates[h_imax].score == 0)
      h_imax = 0;

    //report result
    *selected_cand = imax;
    *selected_h_cand = h_imax;
  }

  bool OptimumSelector::select_refinement(Element* element, int quad_order, Solution* rsln, ElementToRefine& refinement) {
    //make an uniform order in a case of a triangle
    int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
    if (element->is_triangle()) {
      assert_msg(order_v == 0, "Element %d is a triangle but order_v (%d) is not zero", element->id, order_v);
      order_v = order_h;
      quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v); //in a case of a triangle, order_v is zero. Set it to order_h in order to simplify the routines.
    }

    //check validity
    assert_msg(std::max(order_h, order_v) <= H2DRS_MAX_ORDER, "Given order (%d, %d) exceedes the maximum supported order %d.", order_h, order_v, H2DRS_MAX_ORDER);

    //set shapeset mode
    if (element->is_triangle())
      shapeset->set_mode(H2D_MODE_TRIANGLE);
    else
      shapeset->set_mode(H2D_MODE_QUAD);

    //set orders
    set_current_order_range(element);

    //build candidates
    int inx_cand, inx_h_cand;
    create_candidates(element, quad_order
      , H2D_MAKE_QUAD_ORDER(current_max_order, current_max_order)
      , H2D_MAKE_QUAD_ORDER(current_max_order, current_max_order));
    if (candidates.size() > 1) { //there are candidates to choose from
      // evaluate candidates (sum partial projection errors, calculate dofs)
      double avg_error, dev_error;
      evaluate_candidates(element, rsln, &avg_error, &dev_error);

      //select candidate
      select_best_candidate(element, avg_error, dev_error, &inx_cand, &inx_h_cand);
    }
    else { //there is not candidate to choose from, select the original candidate
      inx_cand = 0;
      inx_h_cand = 0;
    }

    //copy result to output
    Cand& cand = candidates[inx_cand];
    Cand& cand_h = candidates[inx_h_cand];
    refinement.split = cand.split;
    ElementToRefine::copy_orders(refinement.p, cand.p);
    if (candidates[inx_h_cand].split == H2D_REFINEMENT_H) { //inx_h_cand points to a candidate which is a H-candidate: copy orders
      ElementToRefine::copy_orders(refinement.q, cand_h.p);
    }
    else { //the index is not H-candidate because not candidate was generate: so, fake orders
      int h_cand_orders[H2D_MAX_ELEMENT_SONS] = { cand_h.p[0], cand_h.p[0], cand_h.p[0], cand_h.p[0] };
      ElementToRefine::copy_orders(refinement.q, h_cand_orders);
    }

    //modify orders in a case of a triangle such that order_v is zero
    if (element->is_triangle())
    {
      for(int i = 0; i < H2D_MAX_ELEMENT_SONS; i++) {
        assert_msg(H2D_GET_V_ORDER(refinement.p[i]) == 0 || H2D_GET_H_ORDER(refinement.p[i]) == H2D_GET_V_ORDER(refinement.p[i]), "Triangle processed but the resulting order (%d, %d) of son %d is not uniform", H2D_GET_H_ORDER(refinement.p[i]), H2D_GET_V_ORDER(refinement.p[i]), i);
        refinement.p[i] = H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(refinement.p[i]), 0);
        assert_msg(H2D_GET_V_ORDER(refinement.q[i]) == 0 || H2D_GET_H_ORDER(refinement.q[i]) == H2D_GET_V_ORDER(refinement.q[i]), "Triangle processed but the resulting q-order (%d, %d) of son %d is not uniform", H2D_GET_H_ORDER(refinement.q[i]), H2D_GET_V_ORDER(refinement.q[i]), i);
        refinement.q[i] = H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(refinement.q[i]), 0);
      }
    }

    if (inx_cand == 0)
      return false;
    else
      return true;
  }

  void OptimumSelector::update_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders) {
    assert_msg(refinement != H2D_REFINEMENT_P, "P-candidate not supported for updating shared orders");
    const int num_sons = get_refin_sons(refinement);
    if (suggested_quad_orders != NULL) {
      for(int i = 0; i < num_sons; i++)
        tgt_quad_orders[i] = suggested_quad_orders[i];
    }
    else {
      //calculate new quad_orders
      int quad_order = orig_quad_order;
      if (cand_list != H2D_H_ISO && cand_list != H2D_H_ANISO) { //H_ISO and H_ANISO has to keep given order
        int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
        switch(refinement) {
          case H2D_REFINEMENT_H:
            order_h = std::max(1, (order_h+1)/2);
            order_v = std::max(1, (order_v+1)/2);
            break;
          case H2D_REFINEMENT_ANISO_H:
            order_v = std::max(1, 2*(order_v+1)/3);
            break;
          case H2D_REFINEMENT_ANISO_V:
            order_h = std::max(1, 2*(order_h+1)/3);
            break;
        }
        if (element->is_triangle())
          quad_order = H2D_MAKE_QUAD_ORDER(order_h, 0);
        else
          quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
      }
      for(int i = 0; i < num_sons; i++)
        tgt_quad_orders[i] = quad_order;
    }
#ifdef _DEBUG
    for(int i = num_sons; i < H2D_MAX_ELEMENT_SONS; i++)
      tgt_quad_orders[i] = 0;
#endif
  }
}
