//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVPressureVariable.h"
#include "ADFunctorInterface.h"

class InputParameters;

/**
 * A variable class for pressure that handle peridic pressure "shift" as
 * extrapolated or Dirichlet boundary faces. The downwind is flagged as an extrapolated boundary
 * face while the upwind face is flagged as a Dirichlet face. 
 * 
 * FIXME: The upwind Dirichlet value is computed
 * using the downwind extrapolated pressure value via linear interpolation using face info's weighting factor _gc
 */
class PeriodicPressureVariable : public INSFVPressureVariable, public ADFunctorInterface
{
public:
  PeriodicPressureVariable(const InputParameters & params);

  static InputParameters validParams();

  bool isExtrapolatedBoundaryFace(const FaceInfo & fi,
                                  const Elem * elem,
                                  const Moose::StateArg & time) const override;

  VectorValue<ADReal> adGradSln(const FaceInfo & fi,
                                const Moose::StateArg & state,
                                bool correct_skewness = false) const override;
                                  
  void initialSetup() override;

protected:
  bool isDirichletBoundaryFace(const FaceInfo & fi,
                               const Elem * elem,
                               const Moose::StateArg & time) const override;

  ADReal getDirichletBoundaryFaceValue(const FaceInfo & fi,
                                       const Elem * elem,
                                       const Moose::StateArg & time) const override;


  /**
   * Checks to see whether the provided element is upwind of the provided face
   * @param fi the face
   * @return a pair in which the first member is whether the element is upwind of the face and the
   * second member is the evaluated face superficial velocity
   */
  bool
  isPeriodic(const FaceInfo & fi) const;

  std::pair<unsigned int, bool>
  periodicPair(const FaceInfo & fi) const;

  Real
  periodicShift(const Elem * me, const unsigned int my_side) const;

  /*
  /// The x-component of velocity
  const Moose::Functor<ADReal> * _u;
  /// The y-component of velocity
  const Moose::Functor<ADReal> * _v;
  /// The z-component of velocity
  const Moose::Functor<ADReal> * _w;
  /// The porosity
  const Moose::Functor<ADReal> * _eps;
  /// The density
  const Moose::Functor<ADReal> * _rho;
  */

  /// The names of the sidesets which will have associated form loss coefficients
  std::vector<BoundaryName> _periodic_sidesets;

  /// The boundary IDs corresponding to the form loss sidesets
  const std::vector<BoundaryID> _periodic_sideset_ids;

  /// Pressure head lost from primary to secondary, one per pair
  const std::vector<Real> _periodic_heads;

};
