// Copyright (c) 2018-2020  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_POLYHEDRON_H
#define CGAL_DRAW_POLYHEDRON_H

#include <CGAL/license/Polyhedron.h>
#include <CGAL/Graphic_storage.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

namespace CGAL
{

#define CGAL_POLY_TYPE CGAL::Polyhedron_3 \
  <PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>

// Specialization of add_in_graphic_storage function.
template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc,
         typename BufferType=float,
         class DrawingFunctor>
void add_in_graphic_storage(const CGAL_POLY_TYPE& apoly,
                           CGAL::Graphic_storage<BufferType> &graphic_storage,
                           const DrawingFunctor &drawing_functor)
{ add_in_graphic_storage_for_fg(apoly, graphic_storage, drawing_functor); }

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc,
         typename BufferType=float>
void add_in_graphic_storage(const CGAL_POLY_TYPE& apoly,
                           CGAL::Graphic_storage<BufferType> &graphic_storage)
{ add_in_graphic_storage_for_fg(apoly, graphic_storage); }

// Specialization of draw function: require Qt and the CGAL basic viewer.
#ifdef CGAL_USE_BASIC_VIEWER

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc,
         typename BufferType=float>
void draw(const CGAL_POLY_TYPE& apoly,
          const char* title="Polyhedron Basic Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  add_in_graphic_storage_for_fg(apoly, buffer);
  draw_graphic_storage(buffer, title);
}

template<class PolyhedronTraits_3,
         class PolyhedronItems_3,
         template < class T, class I, class A>
         class T_HDS,
         class Alloc,
         typename BufferType=float,
         class DrawingFunctor>
void draw(const CGAL_POLY_TYPE& apoly,
          const DrawingFunctor &drawing_functor,
          const char* title="Polyhedron Basic Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  add_in_graphic_storage_for_fg(apoly, buffer, drawing_functor);
  draw_graphic_storage(buffer, title);
}
#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_POLY_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_POLYHEDRON_H
