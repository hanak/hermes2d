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
#include "forms.h"
#include "refmap.h"
#include "integrals_hcurl.h"
#include "adapt.h"
#include "hcurl_adapt.h"

using namespace std;

#ifdef COMPLEX

HcurlAdapt::HcurlAdapt(const Tuple<Space*>& spaces) : Adapt(spaces) {
  for (int i = 0; i < num_comps; i++)
    for (int j = 0; j < num_comps; j++) {
      if (i == j) {
        form[i][j] = hcurl_form<double, scalar>;
        ord[i][j]  = hcurl_form<Ord, Ord>;
      }
    }
}

void HcurlAdapt::prepare_eval_error_value(const int gip_inx, const Func<scalar>& err_sln, const Func<scalar>& rsln) {
  err_sln.val0[gip_inx] = err_sln.val0[gip_inx] - rsln.val0[gip_inx];
  err_sln.val1[gip_inx] = err_sln.val1[gip_inx] - rsln.val1[gip_inx];
  err_sln.curl[gip_inx] = err_sln.curl[gip_inx] - rsln.curl[gip_inx];
}

#endif