//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FaceInfo.h"
#include "MooseTypes.h"

#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/remote_elem.h"

FaceInfo::FaceInfo(const ElemInfo * elem_info,
                   unsigned int side,
                   const dof_id_type id,
                   libMesh::ElemSideBuilder & side_builder)
  : _elem_info(elem_info),
    _neighbor_info(nullptr),
    _id(id),
    _processor_id(_elem_info->elem()->processor_id()),
    _elem_side_id(side),
    _neighbor_side_id(libMesh::invalid_uint),
    _gc(0.5)
{
  // Compute face-related quantities
  auto & face = side_builder(*_elem_info->elem(), side);
  _face_area = face.volume();
  _face_centroid = face.vertex_average();
  _normal = _elem_info->elem()->side_vertex_average_normal(side);
}

void
FaceInfo::computeInternalCoefficients(const ElemInfo * const neighbor_info)
{
  mooseAssert(neighbor_info,
              "We need a neighbor if we want to compute interpolation coefficients!");
  _neighbor_info = neighbor_info;
  _neighbor_side_id = _neighbor_info->elem()->which_neighbor_am_i(_elem_info->elem());
  computeCoefficientsFromNeighborCentroid(_neighbor_info->centroid());

  // Setup quantities used for the approximation of the spatial derivatives
  _d_cn = _neighbor_info->centroid() - _elem_info->centroid();
  _d_cn_mag = _d_cn.norm();
  _e_cn = _d_cn / _d_cn_mag;

  Point r_intersection =
      _elem_info->centroid() +
      (((_face_centroid - _elem_info->centroid()) * _normal) / (_e_cn * _normal)) * _e_cn;

  // For interpolation coefficients
  _gc = (_neighbor_info->centroid() - r_intersection).norm() / _d_cn_mag;
}

void
FaceInfo::computeBoundaryCoefficients()
{
  mooseAssert(!_neighbor_info, "This functions shall only be called on a boundary!");

  // Setup quantities used for the approximation of the spatial derivatives
  _d_cn = _face_centroid - _elem_info->centroid();
  _d_cn_mag = _d_cn.norm();
  _e_cn = _d_cn / _d_cn_mag;

  // For interpolation coefficients
  _gc = 1.0;
}

Point
FaceInfo::skewnessCorrectionVector() const
{
  const Point r_intersection =
      _elem_info->centroid() +
      (((_face_centroid - _elem_info->centroid()) * _normal) / (_e_cn * _normal)) * _e_cn;

  return _face_centroid - r_intersection;
}

void
FaceInfo::computePeriodicCoefficients(const ElemInfo * const neighbor_info,
                                      const Point & translated_centroid,
                                      const unsigned int neighbor_side)
{
  mooseAssert(neighbor_info, "We need a neighbor to compute periodic coefficients");
  mooseAssert(neighbor_side < neighbor_info->elem()->n_sides(), "Invalid neighbor side");
  mooseAssert(neighbor_info->elem() != _elem_info->elem(),
              "A single element spanning the periodic direction is not supported");

  _neighbor_info = neighbor_info;
  _neighbor_side_id = neighbor_side;
  _periodic = true;
  _periodic_neighbor_centroid = translated_centroid;
  _periodic_offset = translated_centroid - neighbor_info->centroid();

  // The image must lie on the outward side of the face; otherwise the translation sign is wrong
  mooseAssert((translated_centroid - _face_centroid) * _normal > 0,
              "Periodic neighbor image on the wrong side of the face (translation direction?)");

  computeCoefficientsFromNeighborCentroid(translated_centroid);
}

void
FaceInfo::computeCoefficientsFromNeighborCentroid(const Point & xN)
{
  const Point & xC = _elem_info->centroid();

  _d_cn = xN - xC;
  _d_cn_mag = _d_cn.norm();
  mooseAssert(_d_cn_mag > 0, "Coincident centroids");
  _e_cn = _d_cn / _d_cn_mag;

  const Real denom = _e_cn * _normal;
  mooseAssert(!MooseUtils::absoluteFuzzyEqual(denom, 0.), "Cell-to-cell line is parallel to the face");

  // intersection of the C-N line with the face plane
  const Point r_intersection = xC + (((_face_centroid - xC) * _normal) / denom) * _e_cn;

  // phi_f = gc * phi_C + (1 - gc) * phi_N
  _gc = (xN - r_intersection).norm() / _d_cn_mag;

  _skewness_correction_vector = _face_centroid - r_intersection;
}