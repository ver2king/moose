//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "PeriodicPressureJump.h"

#include "FVUtils.h"
#include "MooseEnum.h"

registerMooseObject("NavierStokesApp", PeriodicPressureJump);

InputParameters
PeriodicPressureJump::validParams()
{
  auto params = PressureJumpModel::validParams();
  params.addClassDescription(
      "Provides pressure jumps for periodic BCs to PorousRhieChowMassFlux.");
  params.addParam<std::vector<Real>>(
        "pressure_shifts", {}, "Pressure shift through the periodic BCs. ");
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC |
          Moose::RelationshipManagerType::COUPLING,
      [](const InputParameters & /*obj_params*/, InputParameters & rm_params)
      { rm_params.set<unsigned short>("layers") = 3; });
      return params;
}

PeriodicPressureJump::PeriodicPressureJump(const InputParameters & params)
  : PressureJumpModel(params),
    _periodic_heads(getParam<std::vector<Real>>("pressure_shifts"))
{
  if (_boundary_ids.size() != 2 * _periodic_heads.size())
    paramError("pressure_shifts", "Provide one shift per (primary, secondary) boundary pair. ");

  mooseInfo("Dimensions of boundary IDs = ", _boundary_ids.size());
  mooseInfo("Dimensions of pressure shifts = ", _periodic_heads.size());

  for (int i = 0 ; i < _periodic_heads.size(); ++i)
  { 
    auto ids_it = std::next(_boundary_ids.begin(), 2 * i);
    auto next_ids_it = std::next(_boundary_ids.begin(), 2 * i + 1);
    // owner on primary
    _jump_by_owner_id[*ids_it] = -_periodic_heads[i]; 
    // owner on secondary 
    _jump_by_owner_id[*next_ids_it] = _periodic_heads[i]; 
  }
}

Real
PeriodicPressureJump::computePressureJump(const FaceInfo & fi,
                                          const Real /*face_mass_flux*/) const
{

  mooseInfo("Compute pressure jump (periodic BC) ... ");

  //mooseAssert(appliesTo(fi), "The pressure jump model must apply to the supplied face.");
  mooseAssert(fi.isPeriodic(), "A periodic pressure jump requires a periodic face.");

  const bool elem_is_owner = elemIsOwner(fi);
  const Elem * const owner = elem_is_owner ? fi.elemPtr() : fi.neighborPtr();
  const auto owner_side = elem_is_owner ? fi.elemSideID() : fi.neighborSideID();

  Real jump = 0.0;
  unsigned int n_matches = 0;
  for (const auto id : _moose_mesh.getBoundaryIDs(owner, owner_side))
    if (const auto it = _jump_by_owner_id.find(id); it != _jump_by_owner_id.end())
    {
      jump += it->second;
      ++n_matches;
    }

  if (n_matches != 1)
    mooseError("Periodic face ", fi.id(), " matched ", n_matches, " periodic sidesets.");
  
  mooseInfo("Total pressure jump (periodic BCs) = ", jump);
  return jump;
}