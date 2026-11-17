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
  params.addParam<std::vector<BoundaryName>>(
        "periodic_sidesets", {}, "Sidesets over which form the periodic BCs.");
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
    _periodic_sidesets(getParam<std::vector<BoundaryName>>("periodic_sidesets")),
    _periodic_sideset_ids(this->_mesh.getBoundaryIDs(_periodic_sidesets)),
    _periodic_heads(getParam<std::vector<Real>>("pressure_shifts"))
{
  mooseInfo("Dimension of periodic sideset IDs = ", _periodic_sideset_ids.size());
  for (auto id : _periodic_sideset_ids)
  {
    std::cerr << "Periodic sideset id = " <<  id << " is retrieved. " << std::endl;
  };
}

bool
PeriodicPressureJump::isPeriodic(const FaceInfo & fi) const
{   
  // We must ensure any sidesets from the input parameter belong to the face info 
  // so that isPeridic return true. Otherwise return false.
  if (!fi.neighborPtr())
    return false; // this is a boundary that is not periodic

  const auto & binfo = _mesh.getMesh().get_boundary_info();
  for (const auto id : _periodic_sideset_ids)
    if (binfo.has_boundary_id(&fi.elem(), fi.elemSideID(), id) ||
        binfo.has_boundary_id(fi.neighborPtr(), fi.neighborSideID(), id))
      return true;
  return false;
}

Real
PeriodicPressureJump::computePressureJump(const FaceInfo & fi,
                                        const Real face_mass_flux) const
{   
  mooseWarning("PeriodicPressureJump does not use face_mass_flux. ");

  const auto time = Moose::currentState();
  Real pressure_jump = 0.;
  return pressure_jump;
}