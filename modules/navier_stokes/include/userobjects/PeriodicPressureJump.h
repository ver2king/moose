//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "PressureJumpModel.h"

/**
 * Computes the pressure gradient throughout a face for periodic BCs, via FaceInfo
 */
class PeriodicPressureJump : public PressureJumpModel
{
public:
  /// Return the input parameters for this pressure-jump model.
  static InputParameters validParams();

  /**
   * Construct the periodic pressure jump model.
   * @param params The input parameters for the model
   */
  PeriodicPressureJump(const InputParameters & params);

  /**
   * Compute pressure throughout FaceInfo fi
   * FIXME: Should determine the notations for "from what elem to what elem" ?
   * @param fi The internal face carrying the pressure jump
   * @param face_mass_flux The signed mass flux through the face
   */
  Real computePressureJump(const FaceInfo & fi, Real face_mass_flux) const override;

private:

  /**
   * Checks to see whether the provided element is upwind of the provided face
   * @param fi the face
   * @return true if the boundary is periodic, false otherwise 
   */
  bool
  isPeriodic(const FaceInfo & fi) const;

private:
  /// The names of the sidesets which will have associated form loss coefficients
  std::vector<BoundaryName> _periodic_sidesets;

  /// The boundary IDs corresponding to the form loss sidesets
  const std::vector<BoundaryID> _periodic_sideset_ids;

  /// Pressure head lost from primary to secondary, one per pair
  const std::vector<Real> _periodic_heads;
};
