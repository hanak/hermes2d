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

#include "common.h"
#include "limit_order.h"
#include "solution.h"
#include "linsystem.h"
#include "refmap.h"
#include "quad_all.h"
#include "matrix.h"
#include "traverse.h"
#include "norm.h"
#include "element_to_refine.h"
#include "ref_selectors/selector.h"
#include "adapt.h"

using namespace std;

Adapt::Adapt(const Tuple<Space*>& spaces)
  : num_comps(spaces.size()) {
  // check validity
  error_if(num_comps <= 0, "to few components (%d), only %d supported", num_comps, H2D_MAX_COMPONENTS);
  error_if(num_comps >= H2D_MAX_COMPONENTS, "to many components (%d), only %d supported", num_comps, H2D_MAX_COMPONENTS);

  //initialize spaces
  for(int i = 0; i < num_comps; i++)
    this->spaces[i] = spaces[i];

  // reset values
  have_errors = false;
  memset(errors, 0, sizeof(errors));
  memset(form, 0, sizeof(form));
  memset(ord, 0, sizeof(ord));
}

Adapt::~Adapt()
{
  for (int i = 0; i < num_comps; i++)
    if (errors[i] != NULL)
      delete [] errors[i];
}

//// adapt /////////////////////////////////////////////////////////////////////////////////////////

bool Adapt::adapt(RefinementSelectors::Selector* refinement_selector, double thr, int strat, int regularize, bool same_orders, double to_be_processed)
{
  error_if(!have_errors, "element errors have to be calculated first, call calc_error().");
  error_if(refinement_selector == NULL, "selector not provided");

  //get meshes
  int max_id = -1;
  Mesh* meshes[H2D_MAX_COMPONENTS];
  for (int j = 0; j < num_comps; j++) {
    meshes[j] = spaces[j]->get_mesh();
    rsln[j]->set_quad_2d(&g_quad_2d_std);
    rsln[j]->enable_transform(false);
    if (meshes[j]->get_max_element_id() > max_id)
      max_id = meshes[j]->get_max_element_id();
  }

  //reset element refinement info
  AUTOLA2_OR(int, idx, max_id + 1, num_comps + 1);
  for(int j = 0; j < max_id; j++)
    for(int l = 0; l < num_comps; l++)
      idx[j][l] = -1; // element not refined

  double err0 = 1000.0;
  double processed_error = 0.0;

  vector<ElementToRefine> elem_inx_to_proc; //list of indices of elements that are going to be processed
  elem_inx_to_proc.reserve(num_act_elems);

  //adaptivity loop
  double error_threshod = -1; //an error threshold that breaks the adaptivity loop in a case of strategy 1
  int num_exam_elem = 0; //a number of examined elements
  int num_ignored_elem = 0; //a number of ignored elements
  int num_not_changed = 0; //a number of element that were not changed
  int num_priority_elem = 0; //a number of elements that were processed using priority queue

  double error_threshold = 0; //error threshold for stategy 1
  bool first_regular_element = true; //true if first regular element was not processed yet
  int inx_regular_element = 0;
  while (inx_regular_element < num_act_elems || !priority_queue.empty())
  {
    int id, comp, inx_element;

    //get element identification
    if (priority_queue.empty()) {
      id = regular_queue[inx_regular_element].id;
      comp = regular_queue[inx_regular_element].comp;
      inx_element = inx_regular_element;
      inx_regular_element++;
    }
    else {
      id = priority_queue.front().id;
      comp = priority_queue.front().comp;
      inx_element = -1;
      priority_queue.pop();
      num_priority_elem++;
    }
    num_exam_elem++;

    //get info linked with the element
    double err = errors[comp][id];
    Mesh* mesh = meshes[comp];
    Element* e = mesh->get_element(id);

    if (!should_ignore_element(inx_element, mesh, e)) {
      //check if adaptivity loop should end
      if (inx_element >= 0) {
        //prepare error threshold for strategy 1
        if (first_regular_element) {
          error_threshod = thr * err;
          first_regular_element = false;
        }

        // first refinement strategy:
        // refine elements until prescribed amount of error is processed
        // if more elements have similar error refine all to keep the mesh symmetric
        if ((strat == 0) && (processed_error > sqrt(thr) * total_err) && fabs((err - err0)/err0) > 1e-3) break;

        // second refinement strategy:
        // refine all elements whose error is bigger than some portion of maximal error
        if ((strat == 1) && (err < error_threshod)) break;

        if ((strat == 2) && (err < thr)) break;

        if ((strat == 3) &&
          ( (err < error_threshod) ||
          ( processed_error > 1.5 * to_be_processed )) ) break;        
      }

      // get refinement suggestion
      ElementToRefine elem_ref(id, comp);
      int current = spaces[comp]->get_element_order(id);
      bool refined = refinement_selector->select_refinement(e, current, rsln[comp], elem_ref);

      //add to a list of elements that are going to be refined
      if (refined && can_adapt_element(mesh, e, elem_ref.split, elem_ref.p, elem_ref.q) ) {
        idx[id][comp] = (int)elem_inx_to_proc.size();
        elem_inx_to_proc.push_back(elem_ref);
        err0 = err;
        processed_error += err;
      }
      else {
        debug_log("Element (id:%d, comp:%d) not changed", e->id, comp);
        num_not_changed++;
      }
    }
    else {
      num_ignored_elem++;
    }
  }

  verbose("Examined elements: %d", num_exam_elem);
  verbose(" Elements taken from priority queue: %d", num_priority_elem);
  verbose(" Ignored elements: %d", num_ignored_elem);
  verbose(" Not changed elements: %d", num_not_changed);
  verbose(" Elements to process: %d", elem_inx_to_proc.size());
  bool done = false;
  if (num_exam_elem == 0)
    done = true;
  else if (elem_inx_to_proc.empty())
  {
    warn("None of the elements selected for refinement could be refined. Adaptivity step not successful, returning 'true'.");
    done = true;
  }

  //fix refinement if multimesh is used
  fix_shared_mesh_refinements(meshes, num_comps, elem_inx_to_proc, idx, refinement_selector);

  //apply refinements
  apply_refinements(meshes, elem_inx_to_proc);

  //homogenize orders
  if (same_orders)
    homogenize_shared_mesh_orders(meshes);

  // mesh regularization
  if (regularize >= 0)
  {
    if (regularize == 0)
    {
      regularize = 1;
      warn("Total mesh regularization is not supported in adaptivity. 1-irregular mesh is used instead.");
    }
    for (int i = 0; i < num_comps; i++)
    {
      int* parents;
      parents = meshes[i]->regularize(regularize);
      spaces[i]->distribute_orders(meshes[i], parents);
      delete [] parents;
    }
  }

  for (int j = 0; j < num_comps; j++)
    rsln[j]->enable_transform(true);

  verbose("Refined elements: %d", elem_inx_to_proc.size());

  have_errors = false;
  if (strat == 2 && done == true)
    have_errors = true; // space without changes

  return done;
}

void Adapt::fix_shared_mesh_refinements(Mesh** meshes, const int num_comps, std::vector<ElementToRefine>& elems_to_refine, AutoLocalArray2<int>& idx, RefinementSelectors::Selector* refinement_selector) {
  int num_elem_to_proc = elems_to_refine.size();
  for(int inx = 0; inx < num_elem_to_proc; inx++) {
    ElementToRefine& elem_ref = elems_to_refine[inx];
    int current_quad_order = spaces[elem_ref.comp]->get_element_order(elem_ref.id);
    Element* current_elem = meshes[elem_ref.comp]->get_element(elem_ref.id);

    //select a refinement used by all components that share a mesh which is about to be refined
    int selected_refinement = elem_ref.split;
    for (int j = 0; j < num_comps; j++)
    {
      if (selected_refinement == H2D_REFINEMENT_H) break; // iso refinement is max what can be recieved
      if (j != elem_ref.comp && meshes[j] == meshes[elem_ref.comp]) { // if a mesh is shared
        int ii = idx[elem_ref.id][j];
        if (ii >= 0) { // and the sample element is about to be refined by another compoment
          const ElementToRefine& elem_ref_ii = elems_to_refine[ii];
          if (elem_ref_ii.split != selected_refinement && elem_ref_ii.split != H2D_REFINEMENT_P) { //select more complicated refinement
            if ((elem_ref_ii.split == H2D_REFINEMENT_ANISO_H || elem_ref_ii.split == H2D_REFINEMENT_ANISO_V) && selected_refinement == H2D_REFINEMENT_P)
              selected_refinement = elem_ref_ii.split;
            else
              selected_refinement = H2D_REFINEMENT_H;
          }
        }
      }
    }

    //fix other refinements according to the selected refinement
    if (selected_refinement != H2D_REFINEMENT_P)
    {
      //get suggested orders for the selected refinement
      const int* suggested_orders = NULL;
      if (selected_refinement == H2D_REFINEMENT_H)
        suggested_orders = elem_ref.q;

      //update orders
      for (int j = 0; j < num_comps; j++) {
        if (j != elem_ref.comp && meshes[j] == meshes[elem_ref.comp]) { // if components share the mesh
          // change currently processed refinement
          if (elem_ref.split != selected_refinement) {
            elem_ref.split = selected_refinement;
            refinement_selector->update_shared_mesh_orders(current_elem, current_quad_order, elem_ref.split, elem_ref.p, suggested_orders);
          }

          // change other refinements
          int ii = idx[elem_ref.id][j];
          if (ii >= 0) {
            ElementToRefine& elem_ref_ii = elems_to_refine[ii];
            if (elem_ref_ii.split != selected_refinement) {
              elem_ref_ii.split = selected_refinement;
              refinement_selector->update_shared_mesh_orders(current_elem, current_quad_order, elem_ref_ii.split, elem_ref_ii.p, suggested_orders);
            }
          }
          else { // element (of the other comp.) not refined at all: assign refinement
            ElementToRefine elem_ref_new(elem_ref.id, j);
            elem_ref_new.split = selected_refinement;
            refinement_selector->update_shared_mesh_orders(current_elem, current_quad_order, elem_ref_new.split, elem_ref_new.p, suggested_orders);
            elems_to_refine.push_back(elem_ref_new);
          }
        }
      }
    }
  }
}

void Adapt::homogenize_shared_mesh_orders(Mesh** meshes) {
  Element* e;
  for (int i = 0; i < num_comps; i++) {
    for_all_active_elements(e, meshes[i]) {
      int current_quad_order = spaces[i]->get_element_order(e->id);
      int current_order_h = get_h_order(current_quad_order), current_order_v = get_v_order(current_quad_order);

      for (int j = 0; j < num_comps; j++)
        if ((j != i) && (meshes[j] == meshes[i])) // components share the mesh
        {
          int quad_order = spaces[j]->get_element_order(e->id);
          current_order_h = std::max(current_order_h, get_h_order(quad_order));
          current_order_v = std::max(current_order_v, get_v_order(quad_order));
        }

      spaces[i]->set_element_order(e->id, make_quad_order(current_order_h, current_order_v));
    }
  }
}

void Adapt::apply_refinements(Mesh** meshes, std::vector<ElementToRefine>& elems_to_refine)
{
  for (vector<ElementToRefine>::const_iterator elem_ref = elems_to_refine.begin(); elem_ref != elems_to_refine.end(); elem_ref++) // go over elements to be refined
  {
    Element* e;
    e = meshes[elem_ref->comp]->get_element(elem_ref->id);

    if (elem_ref->split == H2D_REFINEMENT_P)
      spaces[elem_ref->comp]->set_element_order(elem_ref->id, elem_ref->p[0]);
    else if (elem_ref->split == H2D_REFINEMENT_H) {
      if (e->active)
        meshes[elem_ref->comp]->refine_element(elem_ref->id);
      for (int j = 0; j < 4; j++)
        spaces[elem_ref->comp]->set_element_order(e->sons[j]->id, elem_ref->p[j]);
    }
    else {
      if (e->active)
        meshes[elem_ref->comp]->refine_element(elem_ref->id, elem_ref->split);
      for (int j = 0; j < 2; j++)
        spaces[elem_ref->comp]->set_element_order(e->sons[ (elem_ref->split == 1) ? j : j+2 ]->id, elem_ref->p[j]);
    }
  }
}


///// Unrefinements /////////////////////////////////////////////////////////////////////////////////

void Adapt::unrefine(double thr)
{

  if (!have_errors)
    error("Element errors have to be calculated first, see calc_error().");

  Mesh* mesh[2];
  mesh[0] = spaces[0]->get_mesh();
  mesh[1] = spaces[1]->get_mesh();


  int k = 0;
  if (mesh[0] == mesh[1]) // single mesh
  {
    Element* e;
    for_all_inactive_elements(e, mesh[0])
    {
      bool found = true;
      for (int i = 0; i < 4; i++)
        if (e->sons[i] != NULL && ((!e->sons[i]->active) || (e->sons[i]->is_curved())))
      { found = false;  break; }

      if (found)
      {
        double sum1 = 0.0, sum2 = 0.0;
        int max1 = 0, max2 = 0;
        for (int i = 0; i < 4; i++)
          if (e->sons[i] != NULL)
        {
          sum1 += errors[0][e->sons[i]->id];
          sum2 += errors[1][e->sons[i]->id];
          int oo = spaces[0]->get_element_order(e->sons[i]->id);
          if (oo > max1) max1 = oo;
          oo = spaces[1]->get_element_order(e->sons[i]->id);
          if (oo > max2) max2 = oo;
        }
        if ((sum1 < thr * errors[regular_queue[0].comp][regular_queue[0].id]) &&
             (sum2 < thr * errors[regular_queue[0].comp][regular_queue[0].id]))
        {
          mesh[0]->unrefine_element(e->id);
          mesh[1]->unrefine_element(e->id);
          errors[0][e->id] = sum1;
          errors[1][e->id] = sum2;
          spaces[0]->set_element_order(e->id, max1);
          spaces[1]->set_element_order(e->id, max2);
          k++; // number of unrefined elements
        }
      }
    }
    for_all_active_elements(e, mesh[0])
    {
      for (int i = 0; i < 2; i++)
        if (errors[i][e->id] < thr/4 * errors[regular_queue[0].comp][regular_queue[0].id])
      {
        int oo = H2D_GET_H_ORDER(spaces[i]->get_element_order(e->id));
        spaces[i]->set_element_order(e->id, std::max(oo - 1, 1));
        k++;
      }
    }
  }
  else // multimesh
  {
    for (int m = 0; m < 2; m++)
    {
      Element* e;
      for_all_inactive_elements(e, mesh[m])
      {
        bool found = true;
        for (int i = 0; i < 4; i++)
          if (e->sons[i] != NULL && ((!e->sons[i]->active) || (e->sons[i]->is_curved())))
        { found = false;  break; }

        if (found)
        {
          double sum = 0.0;
          int max = 0;
          for (int i = 0; i < 4; i++)
            if (e->sons[i] != NULL)
          {
            sum += errors[m][e->sons[i]->id];
            int oo = spaces[m]->get_element_order(e->sons[i]->id);
            if (oo > max) max = oo;
          }
          if ((sum < thr * errors[regular_queue[0].comp][regular_queue[0].id]))
          //if ((sum < 0.1 * thr))
          {
            mesh[m]->unrefine_element(e->id);
            errors[m][e->id] = sum;
            spaces[m]->set_element_order(e->id, max);
            k++; // number of unrefined elements
          }
        }
      }
      for_all_active_elements(e, mesh[m])
      {
        if (errors[m][e->id] < thr/4 * errors[regular_queue[0].comp][regular_queue[0].id])
        {
          int oo = H2D_GET_H_ORDER(spaces[m]->get_element_order(e->id));
          spaces[m]->set_element_order(e->id, std::max(oo - 1, 1));
          k++;
        }
      }
    }
  }
  verbose("Unrefined %d elements.", k);
  have_errors = false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt::set_biform(int i, int j, biform_val_t bi_form, biform_ord_t bi_ord)
{
  error_if(i < 0 || i >= num_comps || j < 0 || j >= num_comps, "invalid component number (%d, %d), max. supported components: %d", i, j, H2D_MAX_COMPONENTS);

  form[i][j] = bi_form;
  ord[i][j] = bi_ord;
}

scalar Adapt::eval_error(biform_val_t bi_fn, biform_ord_t bi_ord,
                             MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2,
                             RefMap *rv1,        RefMap *rv2,        RefMap *rrv1,        RefMap *rrv2)
{
  // determine the integration order
  int inc = (rsln1->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ou = init_fn_ord(rsln1->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(rsln2->get_fn_order() + inc);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = bi_ord(1, &fake_wt, ou, ov, fake_e, NULL);
  int order = rrv1->get_inv_ref_order();
  order += o.get_order();
  limit_order(order);

  ou->free_ord(); delete ou;
  ov->free_ord(); delete ov;
  delete fake_e;

  // eval the form
  Quad2D* quad = sln1->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  Geom<double>* e = init_geom_vol(rrv1, order);
  double* jac = rrv1->get_jacobian(order);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * jac[i];

  // function values and values of external functions
  Func<scalar>* err1 = init_fn(sln1, rv1, order);
  Func<scalar>* err2 = init_fn(sln2, rv2, order);
  Func<scalar>* v1 = init_fn(rsln1, rrv1, order);
  Func<scalar>* v2 = init_fn(rsln2, rrv2, order);

  for (int i = 0; i < np; i++) {
    prepare_eval_error_value(i, *err1, *v1);
    prepare_eval_error_value(i, *err2, *v2);
  }

  scalar res = bi_fn(np, jwt, err1, err2, e, NULL);

  e->free(); delete e;
  delete [] jwt;
  err1->free_fn(); delete err1;
  err2->free_fn(); delete err2;
  v1->free_fn(); delete v1;
  v2->free_fn(); delete v2;

  return res;
}


scalar Adapt::eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord,
                            MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2)
{
  // determine the integration order
  int inc = (rsln1->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ou = init_fn_ord(rsln1->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(rsln2->get_fn_order() + inc);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = bi_ord(1, &fake_wt, ou, ov, fake_e, NULL);
  int order = rrv1->get_inv_ref_order();
  order += o.get_order();
  limit_order(order);

  ou->free_ord(); delete ou;
  ov->free_ord(); delete ov;
  delete fake_e;

  // eval the form
  Quad2D* quad = rsln1->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  Geom<double>* e = init_geom_vol(rrv1, order);
  double* jac = rrv1->get_jacobian(order);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * jac[i];

  // function values
  Func<scalar>* v1 = init_fn(rsln1, rrv1, order);
  Func<scalar>* v2 = init_fn(rsln2, rrv2, order);

  scalar res = bi_fn(np, jwt, v1, v2, e, NULL);

  e->free(); delete e;
  delete [] jwt;
  v1->free_fn(); delete v1;
  v2->free_fn(); delete v2;

  return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Adapt::calc_error(Tuple<Solution*> solutions, Tuple<Solution*> ref_solutions) {
  error_if(solutions.size() != ref_solutions.size(), "a number of solutions (%d) and a number of reference solutions (%d) is not the same.", solutions.size(), ref_solutions.size());
  error_if(solutions.size() != num_comps, "wrong number of solutions (%d), expected %d", solutions.size(), num_comps);

  // obtain solutions
  for (int i = 0; i < num_comps; i++) {
    sln[i] = solutions[i];
    sln[i]->set_quad_2d(&g_quad_2d_std);
    rsln[i] = ref_solutions[i];
    rsln[i]->set_quad_2d(&g_quad_2d_std);
  }

  // prepare multi-mesh traversal and error arrays
  AUTOLA_OR(Mesh*, meshes, 2*num_comps);
  Mesh** ref_meshes = meshes + num_comps;
  AUTOLA_OR(Transformable*, tr, 2*num_comps);
  Transformable** ref_tr = tr + num_comps;
  Traverse trav;
  num_act_elems = 0;
  for (int i = 0; i < num_comps; i++)
  {
    meshes[i] = sln[i]->get_mesh();
    meshes[i+num_comps] = rsln[i]->get_mesh();
    tr[i] = sln[i];
    tr[i+num_comps] = rsln[i];

    num_act_elems += sln[i]->get_mesh()->get_num_active_elements();

    int max = meshes[i]->get_max_element_id();
    if (errors[i] != NULL) delete [] errors[i];
    try { errors[i] = new double[max]; }
    catch(bad_alloc&) { error("unable to allocate space for errors of component %d", i); };
    memset(errors[i], 0, sizeof(double) * max);
  }

  //prepare space for norms
  double total_norm = 0.0;
  AUTOLA_OR(double, norms, num_comps);
  memset(norms, 0, num_comps * sizeof(double));
  double total_error = 0.0;

  //calculate error
  Element** ee;
  trav.begin(2*num_comps, meshes, tr);
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    for (int i = 0; i < num_comps; i++)
    {
      RefMap* rmi = sln[i]->get_refmap();
      RefMap* rrmi = rsln[i]->get_refmap();
      for (int j = 0; j < num_comps; j++)
      {
        RefMap* rmj = sln[j]->get_refmap();
        RefMap* rrmj = rsln[j]->get_refmap();
        double e, t;
        if (form[i][j] != NULL)
        {
          #ifndef H2D_COMPLEX
          e = fabs(eval_error(form[i][j], ord[i][j], sln[i], sln[j], rsln[i], rsln[j], rmi, rmj, rrmi, rrmj));
          t = fabs(eval_norm(form[i][j], ord[i][j], rsln[i], rsln[j], rrmi, rrmj));
          #else
          e = std::abs(eval_error(form[i][j], ord[i][j], sln[i], sln[j], rsln[i], rsln[j], rmi, rmj, rrmi, rrmj));
          t = std::abs(eval_norm(form[i][j], ord[i][j], rsln[i], rsln[j], rrmi, rrmj));
          #endif

          norms[i] += t;
          total_norm  += t;
          total_error += e;
          errors[i][ee[i]->id] += e;
        }
      }
    }
  }
  trav.finish();

  //prepare an ordered list of elements according to an error
  fill_regular_queue(meshes, ref_meshes);

  have_errors = true;
  total_err = total_error/* / total_norm*/;
  return sqrt(total_error / total_norm);
}

void Adapt::fill_regular_queue(Mesh** meshes, Mesh** ref_meshes) {
  //prepare space for queue (it is assumed that it will only grow since we can just split)
  regular_queue.clear();
  if (num_act_elems < (int)regular_queue.capacity()) {
    vector<ElementReference> empty_refs;
    regular_queue.swap(empty_refs); //deallocate
    regular_queue.reserve(num_act_elems); //allocate
  }

  //prepare initial fill
  Element* e;
  vector<ElementReference>::iterator elem_info = regular_queue.begin();
  for (int i = 0; i < num_comps; i++)
    for_all_active_elements(e, meshes[i]) {
      regular_queue.push_back(ElementReference(e->id, i));
//       errors[i][e->id] /= norms[i];
// ??? needed or not ???
// when norms of 2 components are very different it can help (microwave heating)
// navier-stokes on different meshes work only without
    }

  //sort
  std::sort(regular_queue.begin(), regular_queue.end(), compare_elems(errors));
}
