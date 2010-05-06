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
  /// Error of a son of a candidate for various order combinations.
  /** If not noted otherwise, the first index is the horizontal order, the second index is the vertical order.
   *   The maximum allowed order is ::H2DRS_MAX_ORDER + 1 (for some Hermes2d internal reason). */
  typedef double SonProjectionError[H2DRS_MAX_ORDER+2][H2DRS_MAX_ORDER+2];

  /// Projection based selector.
  /** Selects a candidate based on an error. The error is calculated for each element of a candidate separatelly and combined together. */
  class H2D_API ProjBasedSelector : public OptimumSelector {
  public: //API
    virtual ~ProjBasedSelector();
  protected:
    ProjBasedSelector(CandList cand_list, double conv_exp, int max_order, Shapeset* shapeset, const Range<int>& vertex_order, const Range<int>& edge_bubble_order); ///< Constructor.

  protected: //error evaluation
#define H2DRS_VALCACHE_INVALID 0 ///< State of value cache: item contains undefined or invalid value.
#define H2DRS_VALCACHE_VALID 1 ///< State of value cache: item contains a valid value.
    /// An item of a value cache.
    template<typename T>
    struct ValueCacheItem {
      bool is_valid() const { return state != H2DRS_VALCACHE_INVALID; }; ///< Returns true, if value is valid.
      void mark(int new_state = H2DRS_VALCACHE_VALID) { state = new_state; }; ///< Sets state of a value cache.
      void set(T new_value) { value = new_value; }; ///< Sets a new value.
      T get() const { return value; }; ///< Returns a value no matter if the item is valid or not.

      ValueCacheItem(const T& value = 0, const int state = H2DRS_VALCACHE_INVALID) : value(value), state(state) {}; ///< Default constructor. By default, it creates a item that contains invalid value.
    private:
      T value; ///< A value stored in the item.
      int state; ///< A state of the image: ::H2DRS_VALCACHE_INVALID or ::H2DRS_VALCACHE_VALID or any other user-defined value.
    };
    typedef double** ProjMatrixCache[H2DRS_MAX_ORDER+2][H2DRS_MAX_ORDER+2]; ///< A projection matrix cache type. Defines a cache of projection matrices for all possible combination of orders.

    ProjMatrixCache proj_matrix_cache[H2D_NUM_MODES]; ///< An array of projection matrices. Used functions are defined through shape_inx. Index to the array is the size. All matrices are square. If record is NULL, the corresponding matrix has to be calculated.
    /// An array of cached RHS values.
    /** The array is allocated in the constructor, size of the array is equal to the maximum index of a shape function + 1.\
     *  It used used strictly during evalution of proj_calc_err_son(). It is kept here in order to avoid frequent reallocating. */
    ValueCacheItem<scalar>* rhs_cache;

    /// Evaluates error of all candidates. Overloaded method OptimumSelector::evaluate_cands_error.
    virtual void evaluate_cands_error(Element* e, Solution* rsln, double* avg_error, double* dev_error);

    /// Calculates projection errors of an elements of candidates for all permitations of orders.
    /** Errors are not normalized. %Range of orders is defined through parameters info_h, info_h, info_aniso. Overloadable.
     *  \param[in] e An element that is being examined by the selector.
     *  \param[in] info_h Information about H-candidates: range of orders, etc.
     *  \param[in] info_p Information about P-candidates: range of orders, etc.
     *  \param[in] info_aniso Information about ANISO-candidates: range of orders, etc.
     *  \param[in] rsln Reference solution.
     *  \param[out] herr Error of elements of H-candidates of various permutation of orders.
     *  \param[out] perr Error of elements of P-candidates of various permutation of orders.
     *  \param[out] anisoerr Error of elements of ANISO-candidates of various permutation of orders. */
    virtual void calc_projection_errors(Element* e, const CandsInfo& info_h, const CandsInfo& info_p, const CandsInfo& info_aniso, Solution* rsln, SonProjectionError herr[4], SonProjectionError perr, SonProjectionError anisoerr[4]);

    /// Calculate projection errors of an element of an candidate.
    /** The element may span over multiple sub-domain. Element of a candidate is defined in the reference domain.
     *  \param[in] mode Mode: ::H2D_MODE_TRIANGLE or ::H2D_MODE_QUAD.
     *  \param[in] gip_points Integration points (location x, location y, weight) in the reference domain.
     *  \param[in] num_gip_points Number of integration point.
     *  \param[in] num_sub Number of subdomains.
     *  \param[in] sub_domains Subdomains (elements of a reference mesh) that occupy the element of a candidate.
     *  \param[in] sub_trfs Transformation from a reference domain of a subdomain to a reference domain of the element of a candidate.
     *  \param[in] sub_rvals Values at integration points for every subdomain (the first index). Contents of this array (the second index) is defined by user through precalc_ref_solution().
     *  \param[in] coefs_mx An array of coefficients that scale df/dx for each subdomain. A coefficient represents effects of a transformation of the subdomain on df/dx.
     *  \param[in] coefs_my An array of coefficients that scale df/dy for each subdomain. A coefficient represents effects of a transformation of the subdomain on df/dy.
     *  \param[in] info Information about candidates: range of orders, etc.
     *  \param[out] errors Calculated errors for a combination of orders. */
    void calc_error_cand_element(const int mode, double3* gip_points, int num_gip_points, const int num_sub, Element** sub_domains, Trf** sub_trfs, scalar*** sub_rvals, double* coefs_mx, double* coefs_my, const CandsInfo& info, SonProjectionError errors);

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

    /// Returns array of arrays of precalculate values of the reference solution at quadrature points of a son of the element for a given son.
    /** This method is intededed to be overloaded in order to provide all necessary values.
     *  The method have to set active element and quadrature inside.
     *  It is assumed that the returned arrays are allocated and deallocated by the child of this class.
     *  Returned array has to be valid during the life-time of the instance.
     *  The method can assume that the an element is refined into 4 elements (sons) in the reference mesh.
     * \param[in] inx_son An index of a son of an element. The index is in a range [0, 3].
     * \param[in] rsln A reference solution. This is source of a mesh.
     * \param[in] element An element of the coarse solution. An element of both the same geometry and the same ID have to be present in the mesh provided by the reference solution.
     * \param[in] intr_gip_order An order of quadrature integration. The number of quadrature points should be retrieved through a quadrature stored in the paremeter rsln.
     * \return A pointer to 2D array: first index is index of a function expansion (e.g., f, df/dx, etc.), second index is index of GIP. */
    virtual scalar** precalc_ref_solution(int inx_son, Solution* rsln, Element* element, int intr_gip_order) = 0;
    virtual double** build_projection_matrix(Shapeset& shapeset, double3* gip_points, int num_gip_points, const int* shape_inx, const int num_shapes) = 0; ///< Builds a projection matrix.
    virtual scalar evaluate_rsh_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, int shape_inx) = 0; ///< Evaluate a single value of the right side for a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.
    virtual double evaluate_error_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) = 0; ///< Evaluate an error of a projection on a sub-element. Provided GIP are defined on a reference domain. Provided transformation will transform form a reference domain of a sub-element to a reference domain of an element.
  };
}

#endif
