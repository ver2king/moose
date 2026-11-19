//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "PressureJumpModel.h"

#include <unordered_map>

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
  /// Pressure head lost from primary to secondary, one per pair
  const std::vector<Real> _periodic_heads;

  ///
  std::unordered_map<BoundaryID, Real> _jump_by_owner_id;
};
