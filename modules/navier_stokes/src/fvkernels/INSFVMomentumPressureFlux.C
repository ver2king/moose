//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumPressureFlux.h"
#include "INSFVVelocityVariable.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", INSFVMomentumPressureFlux);

InputParameters
INSFVMomentumPressureFlux::validParams()
{
  auto params = FVFluxKernel::validParams();
  params += INSFVMomentumResidualObject::validParams();
  params.addClassDescription("Momentum pressure term eps grad_P, as a flux kernel "
                             "using the divergence theoreom, in the  "
                             "incompressible Navier-Stokes momentum equation.");
  params.addRequiredParam<MooseFunctorName>(NS::pressure, "The pressure");
  //params.addParam<std::vector< BoundaryName >>("periodic_boundary_sidesets", (), "Periodic BC sidesets.");
  return params;
}

INSFVMomentumPressureFlux::INSFVMomentumPressureFlux(const InputParameters & params)
  : FVFluxKernel(params), INSFVMomentumResidualObject(*this), _p(getFunctor<ADReal>(NS::pressure))
  //_periodic_boundary_sides(getParam<std::vector<BoundaryName>>("periodic_boundary_sidesets"))
{ 
  if (!dynamic_cast<INSFVVelocityVariable *>(&_var))
    mooseError("INSFVMomentumPressureFlux may only be used with a Navier-Stokes velocity, "
               "of variable type INSFVSuperficialVelocityVariable.");
}

ADReal
INSFVMomentumPressureFlux::computeQpResidual()
{
  ADReal eps_p_interface;
  // Momentum and porosity domains should match
  const auto & face_type = _face_info->faceType(std::make_pair(_var.number(), _var.sys().number()));
  const bool use_elem = (face_type == FaceInfo::VarFaceNeighbors::ELEM) ||
                        (face_type == FaceInfo::VarFaceNeighbors::BOTH);

  const auto * const elem_ptr = use_elem ? &_face_info->elem() : _face_info->neighborPtr();
  const auto & elem = makeElemArg(elem_ptr);
  const auto state = determineState();

  // FIXME: Whether this makes sense to change or otherwise (for periodic BC, pressure grad is 
  // , modified here or otherwise)
  if (onBoundary(*_face_info) && !isPeriodic(*_face_info))
    eps_p_interface = epsilon()(elem, state) * _p(singleSidedFaceArg(), state);
  else if (isPeriodic(*_face_info))
  { 
    // mooseWarning("NOT IMPLEMENTED: Periodic BC requires interpolation of pressure grad flux. ");
    // TODO: Must obtain the periodic element as neighbor here, via MooseMesh for reference
    const auto * neighbor_ptr = _face_info->neighborPtr();
    const auto & neighbor = makeElemArg(neighbor_ptr);

    const auto * elem_ptr = _face_info->elemPtr();
    const auto & elem = makeElemArg(elem_ptr);

    // FIXME: May need to digest where is left/right of the periodic BCs here before interpolation

    if (elem_ptr == nullptr or neighbor_ptr == nullptr)
      mooseError("Periodic BC mush have both element and its neighbor defined. ");

    // FIXME: Maybe modify the iterpolation here ?
    Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                           eps_p_interface,
                           epsilon()(elem, state) * _p(elem, state),
                           epsilon()(neighbor, state) * _p(neighbor, state),
                           *_face_info,
                           true);
  }
  else
  {
    const auto * neighbor_ptr = use_elem ? _face_info->neighborPtr() : &_face_info->elem();
    const auto & neighbor = makeElemArg(neighbor_ptr);

    Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                           eps_p_interface,
                           epsilon()(elem, state) * _p(elem, state),
                           epsilon()(neighbor, state) * _p(neighbor, state),
                           *_face_info,
                           true);
  }

  return eps_p_interface * _normal(_index);
}
