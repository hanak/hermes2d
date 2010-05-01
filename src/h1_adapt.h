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

#ifndef __H2D_H1_ADAPT_H
#define __H2D_H1_ADAPT_H

/// \brief hp-adaptivity module for H1 spaces.
///
/// H1Adapt is a hp-adaptivity module for continuous elements.
/// Given a reference solution, it provides functions to calculate H1 or
/// energy error estimates, acts as a container for the calculated errors.
/// If not specifie by the used, this class uses the most accurate adaptivity
/// selection algorithm which is slow.
///
class H2D_API H1Adapt : public Adapt {
public:
  H1Adapt(const Tuple<Space*>& spaces); ///< Initializes the class.

protected:
  virtual void prepare_eval_error_value(const int gip_inx, const Func<scalar>& err_sln, const Func<scalar>& rsln); ///< Prepare a value for evaluation of error.
};

#endif
