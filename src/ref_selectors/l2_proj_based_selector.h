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

#ifndef __H2D_REFINEMENT_SELECTORS_L2_PROJ_BASED_SELECTOR_H
#define __H2D_REFINEMENT_SELECTORS_L2_PROJ_BASED_SELECTOR_H

#include "proj_based_selector.h"

namespace RefinementSelectors {

  /// A projection-based selector for L2 space. \ingroup g_selectors
  /** This class is designed to be used with the class L2Adapt.
   *  Since an initialization of the class may take a long time,
   *  it is suggested to create the instance outside the adaptivity
   *  loop. */
  class H2D_API L2ProjBasedSelector : public ProjBasedSelector {
  public: //API
    /// Constructor.
    /** \param[in] cand_list A predefined list of candidates.
     *  \param[in] conv_exp A conversion exponent, see evaluate_cands_score().
     *  \param[in] max_order A maximum order which considered. If ::H2DRS_DEFAULT_ORDER, a maximum order supported by the selector is used, see HcurlProjBasedSelector::H2DRS_MAX_L2_ORDER.
     *  \param[in] user_shapeset A shapeset. If NULL, it will use internal instance of the class L2Shapeset. */
    L2ProjBasedSelector(CandList cand_list = H2D_HP_ANISO, double conv_exp = 1.0, int max_order = H2DRS_DEFAULT_ORDER, L2Shapeset* user_shapeset = NULL);
  protected: //overloads
    /// A function expansion of a function f used by this selector.
    enum LocalFuncExpansion {
      H2D_L2FE_VALUE = 0, ///< A function expansion: f.
      H2D_L2FE_NUM = 1 ///< A total considered function expansion.
    };

    static const int H2DRS_MAX_L2_ORDER; ///< A maximum used order in this L2-space selector.

    scalar* precalc_rvals[H2D_MAX_ELEMENT_SONS][H2D_L2FE_NUM]; ///< Array of arrays of precalculates. The first index is an index of a subdomain, the second index is an index of a function expansion (see enum LocalFuncExpansion).

    /// Sets OptimumSelector::current_max_order and OptimumSelector::current_min_order.
    /** The default order range is [1, 8]. If curved, the upper boundary of the range becomes lower.
     *  Overriden function. For details, see OptimumSelector::set_current_order_range().
     *  \todo Replace calculations inside with calculations that uses symbolic constants instead of fixed numbers.
     *  \todo The original implementation uses subtracts 1 in H1 and Hcurl while L2 subtracts 2. Why? */
    virtual void set_current_order_range(Element* element);

    /// Returns an array of values of the reference solution at integration points.
    /**  Overriden function. For details, see ProjBasedSelector::precalc_ref_solution(). */
    virtual scalar** precalc_ref_solution(int inx_son, Solution* rsln, Element* element, int intr_gip_order); ///< Returns array of arrays of precalculate values of a son.

    /// Builds projection matrix using a given set of shapes.
    /**  Overriden function. For details, see ProjBasedSelector::build_projection_matrix(). */
    virtual double** build_projection_matrix(double3* gip_points, int num_gip_points, const int* shape_inx, const int num_shapes); ///< Builds a projection matrix.

    /// Evaluates a value of the right-hande side in a subdomain.
    /**  Overriden function. For details, see ProjBasedSelector::evaluate_rhs_subdomain(). */
    virtual scalar evaluate_rhs_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx); ///> Evaluate a single value of the right side for a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.

    /// Evaluates an squared error of a projection of an element of a candidate onto subdomains.
    /**  Overriden function. For details, see ProjBasedSelector::evaluate_error_squared_subdomain(). */
    virtual double evaluate_error_squared_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj); ///> Evaluate an error of a projection on a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.

  protected: //defaults
    static L2Shapeset default_shapeset; ///< A default shapeset.
  };
}

#endif
