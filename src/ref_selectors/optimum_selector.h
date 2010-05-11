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
#include "../range.h"
#include "order_permutator.h"
#include "selector.h"

#define H2DRS_ASSUMED_MAX_CANDS 512 ///< An estimated maximum number of candidates. Used for purpose of reserving space. \internal \ingroup g_selectors

//TODO: find out why 20 used used, should'n be there 2*(H2DRS_MAX_ORDER+1)
#define H2DRS_INTR_GIP_ORDER 20 ///< An integration order used to integrate while evaluating a candidate. \internal \ingroup g_selectors
#define H2DRS_MAX_ORDER_INC 2 ///< Maximum increase of an order in candidates. \ingroup g_selectors

#define H2DRS_SCORE_DIFF_ZERO 1E-13 ///< A threshold of difference between scores. Anything below this values is considered zero. \internal \ingroup g_selectors

#define H2DRS_ORDER_ANY -1 ///< Any order. Used as a wildcard to indicate that a given order can by any valid order. \internal \ingroup g_selectors

namespace RefinementSelectors {

  /// Predefined list of candidates. \ingroup g_selectors
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

  /// Returns a string representation of a predefined candidate list. \ingroup g_selectors
  /** Used for debugging and output purposes.
   *  \param cand_list A predefined list of candidates.
   *  \return A string representation of the enum value. */
  extern H2D_API const char* get_cand_list_str(const CandList cand_list);

  /// Returns true if a predefined candidate list may contain candidates that are HP. \ingroup g_selectors
  /** \param cand_list A predefined list of candidates.
   *  \return True if a predefined candidate list may contain candidates that are HP. */
  extern H2D_API bool is_hp(const CandList cand_list);

  /// Returns true if a predefined candidate list may contain candidates with an anisotropic change of orders. \ingroup g_selectors
  /** \param cand_list A predefined list of candidates.
   *  \return True if a predefined candidate list may contain candidates with an anisotropic change of orders. */
  extern H2D_API bool is_p_aniso(const CandList cand_list);

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

    /// A shape index.
    /** Any element order higher than both the vertical and the horizontal direction will use a given shape function. */
    struct ShapeInx {
      int order_h; ///< Order in H direction.
      int order_v; ///< Order in V direction.
      int inx; ///< Index of a shape.
      ShapeType type; ///< Shape type.
      ShapeInx(int order_h, int order_v, int inx, ShapeType type) : order_h(order_h), order_v(order_v), inx(inx), type(type) {};
    };

    Shapeset *shapeset; ///< A shapeset used for projections.

    std::vector<ShapeInx> shape_indices[H2D_NUM_MODES]; ///< Shape indices.
    int max_shape_inx[H2D_NUM_MODES]; ///< Maximum of shape indices.
    int next_order_shape[H2D_NUM_MODES][H2DRS_MAX_ORDER+1]; ///< An index of a shape index of the next order.
    bool has_vertex_shape[H2D_NUM_MODES], has_edge_shape[H2D_NUM_MODES], has_bubble_shape[H2D_NUM_MODES]; ///< True if given type is available in the shapeset.

    /// Adds an index (or indices) of a bubble function of a given order if the index was not used yet.
    /** This function allows to add indices of bubble functions that were not yet added on a quadrilateral.
     *  Since a shapeset allows only to retrieve a list of all possible bubbles based on an order of an element,
     *  this functions allows to create back-mapping table: shape index -> the smallest element order that uses it.
     *  The function assumes that all shapes of a lower element order than a given combinations are already added.
     *  Used by build_shape_indices().
     *  \param[in] order_h A horizontal order of an element.
     *  \param[in] order_v A vertical order of an element.
     *  \param[in,out] used_shape_index A vector of used shapes. If a shape index is present in the map, a shape was already added and it will not be added again.
     *  \paarm[in,out] indices A vector of shape indices. The vector is update by the function. */
    void add_bubble_shape_index(int order_h, int order_v, std::map<int, bool>& used_shape_index, std::vector<ShapeInx>& indices);
    void build_shape_indices(const int mode, const Range<int>& vertex_order, const Range<int>& edge_bubble_order); ///< Build shape indices.
    int calc_num_shapes(int mode, int order_h, int order_v, int allowed_type_mask); ///< Calculates number of shapes of up to given order and of allowd types. If order == H2DRS_ORDER_ANY, any order is allowed.

  public:
    OptimumSelector(CandList cand_list, double conv_exp, int max_order, Shapeset* shapeset, const Range<int>& vertex_order, const Range<int>& edge_bubble_order); /// Contructor. Parameters 'vertex_order' and 'edge_bubble_order' are due to the fact that shapesets return valid index even though given shape is invalid.
    virtual ~OptimumSelector() {};
    virtual bool select_refinement(Element* element, int quad_order, Solution* rsln, ElementToRefine& refinement); ///< Selects refinement.
    virtual void generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders); ///< Updates orders of a refinement in another multimesh component which shares a mesh.
  };

  extern H2D_API std::ostream& operator<<(std::ostream& stream, const OptimumSelector::Cand& cand); ///< Flushes contents of a candidate to a stream. Useful for debug print-outs.
}

#endif
