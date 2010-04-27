// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_REFINEMENT_OPTIMUM_SELECTOR_H
#define __H2D_REFINEMENT_OPTIMUM_SELECTOR_H

#include <ostream>
#include "selector.h"

#define H2DRS_ASSUMED_MAX_CANDS 512 ///< A maximum estimated number of candidates. Used for purpose of reserving a space.

#define H2D_NUM_MODES 2 ///< A number of modes.

//TODO: find out why 20 used used, should'n be there 2*(H2DRS_MAX_ORDER+1)
#define H2DRS_INTR_GIP_ORDER 20 ///< Constant GIP order used by during projection to integrate.
#define H2DRS_MAX_ORDER_INC 2 ///< Maximum increase of order in candidates.

#define H2DRS_SCORE_DIFF_ZERO 1E-13 ///< A threshold of difference between scores which is considered zero.

#define H2D_FN_VALUE  0 ///< Index of a function value.
#define H2D_FN_DX     1 ///< Index of df/dx.
#define H2D_FN_DY     2 ///< Index of df/dy.

#define H2DRS_ORDER_ANY -1 ///< Any order.

namespace RefinementSelectors {

  /// Predefined list of candidates.
  enum CandList {
    H2D_P_ISO, ///< P-candidates only. Orders are modified uniformly.
    H2D_P_ANISO, ///< P-candidates only. Orders are modified non-uniformly.
    H2D_H_ISO, ///< H-candidates only. Orders are not modified.
    H2D_H_ANISO, ///< H- and ANISO-candidates only. Orders are not modified.
    H2D_HP_ISO, ///< H- and P-candidates only. Orders are modified uniformly.
    H2D_HP_ANISO_H, ///< H-, ANISO- and P-candidates. Orders are modified uniformly.
    H2D_HP_ANISO_P, ///< H- and P-candidates only. Orders are modified non-uniformly.
    H2D_HP_ANISO ///< H-, ANISO- and P-candidates. Orders are modified non-uniformly.
  };

  extern H2D_API bool is_hp(const CandList cand_list); ///< Returns true if the candidate list contain candidates that are HP.
  extern H2D_API bool is_p_aniso(const CandList cand_list); ///< Returns true if the candidate list contain candidates with anisotropic change of orders.

  template<typename T>
  class Range { ///< Range of values.
    T lower_bound; ///< Lower boundary.
    T upper_bound; ///< Upper boundary.
    bool empty_range; ///< True if range is empty.
  public:
    Range() : empty_range(true) {};
    Range(const T& lower_bound, const T& upper_bound) : empty_range(lower_bound > upper_bound), lower_bound(lower_bound), upper_bound(upper_bound) {};
    bool empty() const { return empty_range; }; ///< Returns true if range is empty.
    const T& lower() const { return lower_bound; }; ///< Returns lower bound.
    const T& upper() const { return upper_bound; }; ///< Returns upper bound.
    bool is_in_closed(const T& value) const { return (value >= lower_bound && value <= upper_bound); }; ///< Returns true if value is inside the closed range.
    bool is_in_open(const T& value) const { return (value > lower_bound && value < upper_bound); }; ///< Returns true if value is inside the open range.

    static Range<T> make_envelope(const Range<T>& a, const Range<T>& b) { ///< Create an envelope which contains both ranges.
      if (a.empty())
        return b;
      else if (b.empty())
        return a;
      else
        return Range(std::min(a.lower(), b.lower()), std::max(a.upper(), b.upper()));
    };
  };

  class H2D_API OrderPermutator { ///< Permutates orders.
  protected:
    int order_h, order_v; ///< Current order
    int start_order_h, start_order_v; ///< Start orders
    int end_order_h, end_order_v; ///< End orders
    bool iso_p; /// True, if increase is iso.
    int* tgt_quad_order; ///< Target quad order to which the order is written automatically. NULL if none.
  public:
    OrderPermutator () {};
    OrderPermutator (int start_quad_order, int end_quad_order, bool iso_p, int* tgt_quad_order = NULL);
    bool next(); ///< Returns false if there is not next combination.
    void reset(); ///< Resets permutator to the start order.
    int get_order_h() const { return order_h; }; ///< Returns horizontal order.
    int get_order_v() const { return order_v; }; ///< Returns vertical order.
    int get_quad_order() const { return H2D_MAKE_QUAD_ORDER(order_h, order_v); }; ///< Returns quad order.
    int get_start_quad_order() const { return H2D_MAKE_QUAD_ORDER(start_order_h, start_order_v); }; ///< Returns start order.
    int get_end_quad_order() const { return H2D_MAKE_QUAD_ORDER(end_order_h, end_order_v); }; ///< Returns start order.
  };

  class H2D_API OptimumSelector : public Selector { ///< Selector that chooses an optimal candidates based on error decrease per a new DOF.
  public: //candidates
    struct Cand { ///< A candidate.
      double error; ///< Error of this candidate.
      int dofs;  ///< Estimated number of DOFs.
      int split; ///< Operation.
      int p[4]; ///< Encoded orders of sons. If V order is zero, V order is equal to U order.
      double score; ///< Score of candidate.

      Cand() {};
      Cand(const int split, const int order_sons[4])
        : split(split) {
          p[0] = order_sons[0];
          p[1] = order_sons[1];
          p[2] = order_sons[2];
          p[3] = order_sons[3];
      };
      Cand(const int split, const int order_son0, const int order_son1 = 0, const int order_son2 = 0, const int order_son3 = 0)
        : split(split) {
          p[0] = order_son0;
          p[1] = order_son1;
          p[2] = order_son2;
          p[3] = order_son3;
      };
      int get_num_sons() const { ///< Returns number of sons.
        switch (split) {
          case H2D_REFINEMENT_H: return 4;
          case H2D_REFINEMENT_P: return 1;
          case H2D_REFINEMENT_ANISO_H:
          case H2D_REFINEMENT_ANISO_V:
            return 2;
          default:
            error("invalid split type %d", split);
            return -1;
            break;
        }
      }

      friend H2D_API std::ostream& operator<<(std::ostream& stream, const Cand& cand);
    };
    const std::vector<Cand>& get_candidates() const { return candidates; }; ///< Returns current candidates.

  protected: //candidates
    struct CandsInfo { ///< Information about candidates.
      bool uniform_orders; ///< True if all sons of all examined candidates have uniform orders
      int min_quad_order; ///< Minimum quad order of all sons of all examined candidates.
      int max_quad_order; ///< Maximum quad order of all sons of all examined candidates.
      CandsInfo() : uniform_orders(false), min_quad_order(-1), max_quad_order(-1) {}; ///< Default constructor. Creates a state that defines no candidates at all.
      bool is_empty() const { return (min_quad_order < 0 || max_quad_order < 0); }; ///< Returns true if there are no examined candidates.
    };

    CandList cand_list; ///< Allowed candidate types.
	  double conv_exp; ///< Convergence power. Modifies difference between DOFs before they are used to calculate the score.
    std::vector<Cand> candidates; ///< A vector of candidates. The first candidate is the original element.
    void update_cands_info(CandsInfo& info_h, CandsInfo& info_p, CandsInfo& info_aniso) const; ///< Updates information (min order, max order) about candidates. Initial information is provided.

    void append_candidates_split(const int start_quad_order, const int last_order, const int split, bool uniform); ///< Creates candidates of a given split-type.

    virtual void create_candidates(Element* e, int quad_order, int max_ha_quad_order, int max_p_quad_order); ///< Initializes the array of candidates and fills it with candidates. Maximum orderds for H- and ANISO-candidates are given by max_ha_quad_order and for P-candidates max_p_quad_order. In the case of ANIO-candidates, the maximum is applied only only to modified orders.
    virtual void evaluate_candidates(Element* e, Solution* rsln, double* avg_error, double* dev_error); ///< Evaluates candidates. Calculates their error and dofs. Calculates average error and sample deviation.
    virtual void select_best_candidate(Element* e, const double avg_error, const double dev_error, int* selected_cand, int* selected_h_cand); ///< Selects the best candidate and the best h-candidate.

    virtual void evaluate_cands_error(Element* e, Solution* rsln, double* avg_error, double* dev_error) = 0; ///< Calculates error of candidates.
    virtual void evaluate_cands_dof(Element* e, Solution* rsln); ///< Calculates DOFs of candidates.

  private:
    static bool compare_cand_score(const Cand& a, const Cand& b); ///< True if score a is greater than score b

  protected: //orders and their range
    int current_max_order; ///< Current maximum order.
    int current_min_order; ///< Current minimum order.

    virtual void set_current_order_range(Element* element) = 0; ///< Sets current maximum and minimum order.

  protected: //shape functions
    enum ShapeType { ///< Shape type.
      H2DST_VERTEX = 0x01, ///< Vertex function.
      H2DST_HORIZ_EDGE = 0x02, ///< Horizontal edge function.
      H2DST_VERT_EDGE = 0x04, ///< Verical edge function.
      H2DST_TRI_EDGE = 0x08, ///< Triangle edge. 
      H2DST_BUBBLE = 0x10 ///< Bubble function.
    };

    struct ShapeInx { ///< A shape index.
      int order_h; ///< Order in H direction. Zero if the shape is just along V-direction.
      int order_v; ///< Order in H direction. Zero if the shape is just along V-direction.
      int inx; ///< Index of a shape.
      ShapeType type; ///< Shape type.
      ShapeInx(int order_h, int order_v, int inx, ShapeType type) : order_h(order_h), order_v(order_v), inx(inx), type(type) {};
    };

    Shapeset *shapeset; ///< A shapeset used for projections.

    std::vector<ShapeInx> shape_indices[H2D_NUM_MODES]; ///< Shape indices.
    int max_shape_inx[H2D_NUM_MODES]; ///< Maximum of shape indices.
    int next_order_shape[H2D_NUM_MODES][H2DRS_MAX_ORDER+1]; ///< An index of a shape index of the next order.
    bool has_vertex_shape[H2D_NUM_MODES], has_edge_shape[H2D_NUM_MODES], has_bubble_shape[H2D_NUM_MODES]; ///< True if given type is available in the shapeset.

    void build_shape_indices(const int mode, const Range<int>& vertex_order, const Range<int>& edge_bubble_order); ///< Build shape indices.
    int calc_num_shapes(int mode, int order_h, int order_v, int allowed_type_mask); ///< Calculates number of shapes of up to given order and of allowd types. If order == H2DRS_ORDER_ANY, any order is allowed.

  public:
    OptimumSelector(CandList cand_list, double conv_exp, int max_order, Shapeset* shapeset, const Range<int>& vertex_order, const Range<int>& edge_bubble_order); /// Contructor. Parameters 'vertex_order' and 'edge_bubble_order' are due to the fact that shapesets return valid index even though given shape is invalid.
    virtual ~OptimumSelector() {};
    virtual bool select_refinement(Element* element, int quad_order, Solution* rsln, ElementToRefine& refinement); ///< Selects refinement.
    virtual void update_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders); ///< Updates orders of a refinement in another multimesh component which shares a mesh.
  };

  extern H2D_API std::ostream& operator<<(std::ostream& stream, const OptimumSelector::Cand& cand); ///< Flushes contents of a candidate to a stream. Useful for debug print-outs.
}

#endif
