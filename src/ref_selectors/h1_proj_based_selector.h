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

#ifndef __HERMES2D_REFINEMENT_SELECTORS_H1_PROJ_BASED_SELECTOR_H
#define __HERMES2D_REFINEMENT_SELECTORS_H1_PROJ_BASED_SELECTOR_H

#include "proj_based_selector.h"

namespace RefinementSelectors {

  class HERMES2D_API H1ProjBasedSelector : public ProjBasedSelector { ///< Selector that does HP-adaptivity in H1 space using projections.
  public: //API
    H1ProjBasedSelector(CandList cand_list = H2D_HP_ANISO, double conv_exp = 1.0, int max_order = H2DRS_DEFAULT_ORDER, H1Shapeset* user_shapeset = NULL);
  protected: //overloads
    scalar* precalc_rvals[H2D_MAX_ELEMENT_SONS][3]; ///< Array of arrays of precalculates values: VALUE, DX, DY.

    virtual void set_current_order_range(Element* element); ///< Sets current maximum and minimum order. If the max_order is H2DRS_DEFAULT_ORDER, in the case of linear elements it uses 9 and in the case of curvilinear elements it depends on iro_cache (how curved they are).

    virtual scalar** precalc_ref_solution(int inx_son, Solution* rsln, Element* element, int intr_gip_order); ///< Returns array of arrays of precalculate values of a son.
    virtual double** build_projection_matrix(Shapeset& shapeset, double3* gip_points, int num_gip_points, const int* shape_inx, const int num_shapes); ///< Builds a projection matrix.
    virtual scalar evaluate_rsh_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx); ///> Evaluate a single value of the right side for a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.
    virtual double evaluate_error_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj); ///> Evaluate an error of a projection on a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.

  protected: //defaults
    static H1Shapeset default_shapeset; ///< Default shapeset.
  };
}

#endif
