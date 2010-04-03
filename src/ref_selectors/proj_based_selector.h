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

#ifndef __H2D_REFINEMENT_PROJ_BASED_SELECTOR_H
#define __H2D_REFINEMENT_PROJ_BASED_SELECTOR_H

#include "../tuple.h"
#include "optimum_selector.h"

namespace RefinementSelectors {
  typedef double SonProjectionError[H2DRS_MAX_ORDER+2][H2DRS_MAX_ORDER+2]; ///< Error of a son of a candidate for various order combinations. The maximum allowed order is H2DRS_MAX_ORDER+1.

  class H2D_API ProjBasedSelector : public OptimumSelector {
  public: //API
    virtual ~ProjBasedSelector();
  protected:
    ProjBasedSelector(AdaptType adapt_type, double conv_exp, int max_order, Shapeset* shapeset, int min_edge_bubble_order);

  protected: //error evaluation
    template<typename T>
    struct ValueCacheItem { ///< An item of a value cache.
      inline bool is_valid() const { return state != 0; }; ///< Returns true, if value is valid.
      inline void mark_invalid() { state = 0; }; ///< Marks a value as invalid.
      inline void mark(int new_state = 1) { state = new_state; }; ///< Sets state of a value cache.
      inline void set(T new_value) { value = new_value; }; ///< Sets a new value.
      inline T get() { return value; }; ///< Returns value.
    private:
      T value; ///< Value.
      int state; ///< State.
    };

    double** proj_matrices[H2DRS_MAX_ORDER+1][H2DRS_MAX_ORDER+1]; ///< An array of projection matrices. Used functions are defined through shape_inx. Index to the array is the size. All matrices are square. If record is NULL, the corresponding matrix has to be calculated.
    ValueCacheItem<scalar>* rhs_cache; ///< An array of RHS values. Valid only during evalution of proj_calc_err_son.

    virtual void evaluate_cands_error(Element* e, Solution* rsln, double* avg_error, double* dev_error); ///< Calculates error of candidates.

    /// \brief Calculate various projection errors for sons of a candidates of given combination of orders. Errors are not normalized. Overloadable.
    /// \param[in] e Element which is being processed.
    /// \param[in] rsln Reference solution.
    /// \param[out] herr Errors of sons of H-candidate. < 0 if there are no H-candidates.
    /// \param[out] anisoerr Errors of sons of ANISO-candidates. < 0 if there are no ANISO-candidates.
    /// \param[out] perr Errros of sons of P-candidates. < 0 if there are no P-candidates.
    virtual void calc_projection_errors(Element* e, const int max_quad_order_h, const int max_quad_order_p, const int max_quad_order_aniso, Solution* rsln, SonProjectionError herr[4], SonProjectionError anisoerr[4], SonProjectionError perr);
    void calc_proj_error_cand_son(const int mode, double3* gip_points, int num_gip_points, const int num_sub, Element** sub_elems, Trf** sub_trfs, scalar*** sub_rvals, double* coefs_mx, double* coefs_my, int max_quad_order, SonProjectionError errors); ///< Calculate projection errors.

  protected: //projection
    struct ElemProj { ///< Element projection parameters.
      int* shape_inxs; ///< Shape indices
      int num_shapes; ///< Number of shape indices.
      scalar* shape_coefs; ///< Coefficients of shape indices of a projection.
      int max_quad_order; ///< Maximum quad order of the projection.
    };
    struct ElemGIP { ///< GIP on a element.
      double3* gip_points; ///< GIP points and weights.
      int num_gip_points; ///< A number of GIP points.
      scalar** rvals; ///< Values of a reference solution at GIP.
    };
    struct ElemSubTrf { ///< Transformation from a coordinate system of a sub-element to coordinate system of an element.
      Trf* trf; ///< Transformation.
      double coef_mx, coef_my; ///< Differentials correction coefficients.
    };

    virtual scalar** precalc_ref_solution(int inx_son, Solution* rsln, Element* element, int intr_gip_order) = 0; ///< Returns array of arrays of precalculate values of a son. First index is index of a value, second index is index of GIP. Before entering, neither son is made active nor order of GIP is set even though both are provied.
    virtual double** build_projection_matrix(Shapeset& shapeset, double3* gip_points, int num_gip_points, const int* shape_inx, const int num_shapes) = 0; ///< Builds a projection matrix.
    virtual scalar evaluate_rsh_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx) = 0; ///> Evaluate a single value of the right side for a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.
    virtual double evaluate_error_sub_element(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) = 0; ///> Evaluate an error of a projection on a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.
  };
}

#endif
