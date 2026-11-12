//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PeriodicPressureVariable.h"
#include "NS.h"
#include "NSFVUtils.h"

registerMooseObject("NavierStokesApp", PeriodicPressureVariable);

InputParameters
PeriodicPressureVariable::validParams()
{
  auto params = INSFVPressureVariable::validParams();
  params += ADFunctorInterface::validParams();
  /*
  params.addRequiredParam<MooseFunctorName>("u", "The x-component of velocity");
  params.addParam<MooseFunctorName>("v", 0, "The y-component of velocity");
  params.addParam<MooseFunctorName>("w", 0, "The z-component of velocity");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "The porosity");
  params.addRequiredParam<MooseFunctorName>(NS::density, "The density");
  */
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

PeriodicPressureVariable::PeriodicPressureVariable(const InputParameters & params)
  : INSFVPressureVariable(params), ADFunctorInterface(this),
    _periodic_sidesets(getParam<std::vector<BoundaryName>>("periodic_sidesets")),
    _periodic_sideset_ids(this->_mesh.getBoundaryIDs(_periodic_sidesets)),
    _periodic_heads(getParam<std::vector<Real>>("pressure_shifts"))
{ 
  mooseInfo("Dimension of periodic sideset IDs = ", _periodic_sideset_ids.size());
  for (auto id : _periodic_sideset_ids)
  {
    std::cerr << "Periodic sideset id = " <<  id << " is retrieved. " << std::endl;
  }
  //mooseInfo("[PeriodicPressureVariable].PeriodicPressureVariable: Variable's kernel. ");
}

void
PeriodicPressureVariable::initialSetup()
{ 
  INSFVPressureVariable::initialSetup();
  
  /// FIXME: Code below is fto debug the variable. 
  // Will be removed after debugging is completed. 
  _console << "[PeriodicPressureVariable] sideset ids:";
  for (const auto id : _periodic_sideset_ids)
    _console << " " << id;
  _console << "  heads: " << _periodic_heads.size() << std::endl;

  unsigned int n_tagged_internal = 0, n_periodic = 0, n_tagged_boundary = 0;
  for (const auto * fi : _mesh.faceInfo())
  {
    const bool tagged = !fi->boundaryIDs().empty();
    if (fi->neighborPtr() && tagged) ++n_tagged_internal;
    if (!fi->neighborPtr() && tagged) ++n_tagged_boundary;
    if (isPeriodic(*fi)) ++n_periodic;
  }
  _console << "[PeriodicPressureVariable] tagged internal faces: " << n_tagged_internal
           << ", tagged boundary faces: " << n_tagged_boundary
           << ", periodic faces: " << n_periodic << std::endl;  
}

bool
PeriodicPressureVariable::isPeriodic(const FaceInfo & fi) const
{   
  // We must ensure any sidesets from the input parameter belong to the face info 
  // so that isPeridic return true. Otherwise return false.
  if (!fi.neighborPtr())
    return false; // no partner: a real boundary, not a seam

  const auto & binfo = _mesh.getMesh().get_boundary_info();
  for (const auto id : _periodic_sideset_ids)
    if (binfo.has_boundary_id(&fi.elem(), fi.elemSideID(), id) ||
        binfo.has_boundary_id(fi.neighborPtr(), fi.neighborSideID(), id))
      return true;
  return false;
}

std::pair<unsigned int, bool>
PeriodicPressureVariable::periodicPair(const FaceInfo & fi) const
{
  // The FaceInfo carries the tags of the side belonging to its owning element,
  // so the tag tells us which sideset fi.elem() sits on.
  for (const auto bd_id : fi.boundaryIDs())
    for (const auto i : index_range(_periodic_sideset_ids))
      if (_periodic_sideset_ids[i] == bd_id)
        // even index = primary, odd = secondary
        return {i / 2, (i % 2 == 0)};

  // Otherwise return error for not having periodic BCs
  mooseError("periodicPair called on a face that does not have periodic BCs. ");
}

Real
PeriodicPressureVariable::periodicShift(const Elem * me, const unsigned int my_side) const
{
  const auto & binfo = _mesh.getMesh().get_boundary_info();
  Real shift = 0.0;
  unsigned int n_matches = 0;
  for (const auto i : index_range(_periodic_sideset_ids))
    if (binfo.has_boundary_id(me, my_side, _periodic_sideset_ids[i]))
    {
      shift += ((i % 2 == 0) ? +1.0 : -1.0) * _periodic_heads[i / 2];
      ++n_matches;
    }
  if (n_matches != 1)
    mooseError("Periodic face side matched ", n_matches, " periodic sidesets; expected 1.");
  return shift;
}

bool
PeriodicPressureVariable::isExtrapolatedBoundaryFace(const FaceInfo & fi,
                                                      const Elem * const elem,
                                                      const Moose::StateArg & time) const
{  
  /*
  if (isDirichletBoundaryFace(fi, elem, time))
    return false;
  if (!isInternalFace(fi))
    // We are neither a Dirichlet nor an internal face
    return true;

  if (isSeparatorBoundary(fi))
    return true;

  if (std::get<0>(NS::isPorosityJumpFace(*_eps, fi, time)))
    // We choose to extrapolate for the element that is downwind
    return !std::get<0>(elemIsUpwind(*elem, fi, time));

  return false;
  */
  return isPeriodic(fi) ? false : INSFVPressureVariable::isExtrapolatedBoundaryFace(fi, elem, time);
  }

bool
PeriodicPressureVariable::isDirichletBoundaryFace(const FaceInfo & fi,
                                                   const Elem * const elem,
                                                   const Moose::StateArg & time) const
{  
  /*
  if (INSFVPressureVariable::isDirichletBoundaryFace(fi, elem, time))
    return true;

  if (isInternalFace(fi))
  {
    if (isSeparatorBoundary(fi))
      return false;

    if (std::get<0>(NS::isPorosityJumpFace(*_eps, fi, time)))
      // We choose to apply the Dirichlet condition for the upwind element
      return std::get<0>(elemIsUpwind(*elem, fi, time));
  }

  return false; 
  */

  if (isPeriodic(fi))
    return elem != nullptr;   // side-aware callers get the jump; two-sided callers interpolate
  return INSFVPressureVariable::isDirichletBoundaryFace(fi, elem, time);
}

ADReal
PeriodicPressureVariable::getDirichletBoundaryFaceValue(const FaceInfo & fi,
                                                         const Elem * const elem,
                                                         const Moose::StateArg & time) const
{ 
  mooseInfo("[PeriodicPressureVariable.getDirichletBoundaryFaceValue]: variable's kernel. ");

  if (INSFVPressureVariable::isDirichletBoundaryFace(fi, elem, time))
    return INSFVPressureVariable::getDirichletBoundaryFaceValue(fi, elem, time);

  // Extract the weighting factor from the corresponding face info to compute pressure shift(s)
  const Elem * const me = elem ? elem : &fi.elem();          // <-- must come first
  mooseAssert(fi.neighborPtr(), "periodic face requires a neighbor. ");

  const bool from_elem = (me == &fi.elem());
  const Elem * const other = from_elem ? fi.neighborPtr() : &fi.elem();
  const unsigned int my_side = from_elem ? fi.elemSideID() : fi.neighborSideID();
  const Real gc = from_elem ? fi.gC() : 1.0 - fi.gC();

  // Accumulate the signed shift over every periodic tag this face carries,
  // in the same spirit as BernoulliPressureVariable's form-factor loop.
  const auto & binfo = _mesh.getMesh().get_boundary_info();
  Real shift = 0.0;
  unsigned int n_matches = 0;
  for (const auto i : index_range(_periodic_sideset_ids))
    if (binfo.has_boundary_id(elem, my_side, _periodic_sideset_ids[i]))
    {
      const bool me_on_primary = (i % 2 == 0);
      // p(secondary) = p(primary) - head: bring the far cell into my frame
      shift += (me_on_primary ? +1.0 : -1.0) * _periodic_heads[i / 2];
      ++n_matches;
    }

  if (n_matches != 1)
    mooseError("Periodic face side matched ", n_matches, " periodic sidesets; expected exactly 1.");

  const ADReal p_me = getElemValue(elem, time);
  const ADReal p_other = getElemValue(other, time) + shift;
  return gc * p_me + (1.0 - gc) * p_other;
}

VectorValue<ADReal>
PeriodicPressureVariable::adGradSln(const FaceInfo & fi,
                                    const Moose::StateArg & time,
                                    const bool correct_skewness) const
{
  if (!isPeriodic(fi))
    return INSFVPressureVariable::adGradSln(fi, time, correct_skewness);

  const Elem * const elem = &fi.elem();
  const Elem * const neighbor = fi.neighborPtr();

  // normal (over-relaxed) difference, with the far cell brought into elem's frame
  const Real shift = periodicShift(elem, fi.elemSideID());
  const ADReal p_e = getElemValue(elem, time);
  const ADReal p_n = getElemValue(neighbor, time) + shift;
  const ADReal normal_grad = (p_n - p_e) / fi.dCNMag();

  // tangential part from the interpolated cell gradients, as the base class does
  const VectorValue<ADReal> grad_e =
      INSFVPressureVariable::adGradSln(elem, time, correct_skewness);
  const VectorValue<ADReal> grad_n =
      INSFVPressureVariable::adGradSln(neighbor, time, correct_skewness);
  VectorValue<ADReal> interp_grad;
  Moose::FV::interpolate(Moose::FV::InterpMethod::Average, interp_grad, grad_e, grad_n, fi, true);

  const auto & e_cn = fi.eCN();
  return interp_grad + (normal_grad - interp_grad * e_cn) * e_cn;
}