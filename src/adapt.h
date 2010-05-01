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

#ifndef __H2D_ADAPT_H
#define __H2D_ADAPT_H

#include "tuple.h"
#include "forms.h"
#include "weakform.h"
#include "integrals_h1.h"
#include "ref_selectors/selector.h"

#define H2D_MAX_COMPONENTS 10 ///< A maximum number of equations. Equals to a mesh-space pairs.

H2D_API_USED_TEMPLATE(Tuple<Space*>);
H2D_API_USED_TEMPLATE(Tuple<Solution*>);

/// \brief hp-adaptivity module.
///
/// H1AdaptHP is a hp-adaptivity module for continuous elements.
/// Given a reference solution, it provides functions to calculate H1 or
/// energy error estimates, acts as a container for the calculated errors.
/// If not specifie by the used, this class uses the most accurate adaptivity
/// selection algorithm which is slow.
///
class H2D_API Adapt
{
protected:
  Adapt(const Tuple<Space*>& spaces); ///< Initializes the class.

public:
  virtual ~Adapt();

  typedef scalar (*biform_val_t) (int n, double *wt, Func<scalar> *u, Func<scalar> *v, Geom<double> *e, ExtData<scalar> *);
  typedef Ord (*biform_ord_t) (int n, double *wt, Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *);

  /// Sets user defined bilinear form to calculate error. Default forms are errors of given norm (on diagonal) are set by default. Use this function only to change it (e.g. energy error).
  void set_biform(int i, int j, biform_val_t bi_form, biform_ord_t bi_ord);

  virtual double calc_error(Tuple<Solution*> solutions, Tuple<Solution*> ref_solutions); /// Calculates error between solution and reference solution.

  /// Refines elements based on results from calc_error() or calc_energy_error().
  bool adapt(RefinementSelectors::Selector* refinement_selector, double thr, int strat = 0,
             int regularize = -1,
             bool same_orders = false, double to_be_processed = 0.0);

  /// Unrefines the elements with the smallest error
  void unrefine(double thr);

  /// Internal. Functions to obtain errors of individual elements.
  struct ElementReference { ///< A reference to a element.
    int id, comp;
    ElementReference() {};
    ElementReference(int id, int comp) : id(id), comp(comp) {};
  };
  double get_element_error(int component, int id) const { error_if(!have_errors, "element errors have to be calculated first, call calc_error()"); return errors[component][id]; } ///< Returns error of an element.
  const std::vector<ElementReference>& get_regular_queue() const { return regular_queue; } ///< Returns regular queue of elements.
  int get_total_active_elements() const { return num_act_elems; } ///< Returns total number of active elements.

  void apply_refinement(const ElementToRefine& elem_ref); ///< Apply refinement.
  virtual void apply_refinements(std::vector<ElementToRefine>& elems_to_refine); ///< Apply refinements.
  const std::vector<ElementToRefine>& get_last_refinements() const { return last_refinements; }; ///< Returns last refinements.

protected: //adaptivity
  int num_act_elems; ///< A total number of active elements across all provided meshes.
  std::queue<ElementReference> priority_queue; ///< A queue of priority elements that are processed before the next element in esort is processed.
  std::vector<ElementReference> regular_queue; ///< A queue of elements which shoulb be examined.
  std::vector<ElementToRefine> last_refinements; ///< A vector of last refinements.

  virtual bool should_ignore_element(const int inx_element, const Mesh* mesh, const Element* element) { return false; }; ///< Returns true, if an element should be ignored for purposes of adaptivity.
  virtual bool can_adapt_element(Mesh* mesh, Element* e, const int split, const int4& p, const int4& q) { return true; }; ///< Returns true, if an element can be adapted using a selected candidate.
  
  void fix_shared_mesh_refinements(Mesh** meshes, const int num_comps, std::vector<ElementToRefine>& elems_to_refine, AutoLocalArray2<int>& idx, RefinementSelectors::Selector* refinement_selector); ///< Fixes refinements of a mesh which is shared among multiple components of a multimesh.
  void homogenize_shared_mesh_orders(Mesh** meshes); ///< Homogenize element orders to the highest order in the mesh if the mesh is shared among components.

protected: // spaces & solutions
  const int num_comps; ///< Number of components.
  Space* spaces[H2D_MAX_COMPONENTS];
  Solution* sln[H2D_MAX_COMPONENTS];
  Solution* rsln[H2D_MAX_COMPONENTS];

protected: // element error arrays
  double* errors[H2D_MAX_COMPONENTS];
  double  norms[H2D_MAX_COMPONENTS];
  bool    have_errors;
  double  total_err;

protected: //forms and error evaluation
  biform_val_t form[H2D_MAX_COMPONENTS][H2D_MAX_COMPONENTS]; ///< Bilinear forms to calculate error
  biform_ord_t ord[H2D_MAX_COMPONENTS][H2D_MAX_COMPONENTS]; ///< Bilinear forms to calculate error

  scalar eval_error(biform_val_t bi_fn, biform_ord_t bi_ord,
                    MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2,
                    RefMap *rv1,        RefMap *rv2,        RefMap *rrv1,        RefMap *rrv2); ///< Evaluate error.

  scalar eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord,
                   MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2); ///< Evalute norm.

  virtual void prepare_eval_error_value(const int gip_inx, const Func<scalar>& err_sln, const Func<scalar>& rsln) = 0; ///< Prepare a value for evaluation of error. The results should be stored to err_sln.


  virtual void fill_regular_queue(Mesh** meshes, Mesh** ref_meshes); ///< Builds a queue of elements to be examined, i.e., inits Adapt::standard_queue. This sorted by error descending. Assumes that Adapt::errors is initialized.
private: 
  class compare_elems { ///< Functor that compares elements. Used to sort elements according to an error using an external array.
  private:
    double** errors;
  public:
    compare_elems(double** errors): errors(errors) {};
    bool operator ()(const ElementReference& e1,const ElementReference& e2) const {
      return errors[e1.comp][e1.id] > errors[e2.comp][e2.id];
    };
  };
};

#endif
