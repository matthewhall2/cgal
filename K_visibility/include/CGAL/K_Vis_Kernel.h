#ifndef MYKERNEL_H
#define MYKERNEL_H

#include <CGAL/Cartesian.h>
#include "K_Vis_point_2.h"
#include "K_VisSegment_2.h"
#include "K_VisConstruct_bbox_2.h"
#include "K_VisConstruct_coord_iterator.h"
#include "K_VisConstruct_point_2.h"

// K_ is the new kernel, and K_Base is the old kernel
template < typename K_, typename K_Base >
class K_Vis_Base
  : public K_Base::template Base<K_>::Type
{
  typedef typename K_Base::template Base<K_>::Type   OldK;
public:
  typedef typename K_                                Kernel;
  typedef typename K_Vis_point_2<OldK>                         Point_2;
  typedef typename K_VisSegment_2<Kernel>               Segment_2;
  typedef typename K_VisConstruct_point_2<Kernel, OldK>       Construct_point_2;
  typedef typename const K_Base::FT*                     Cartesian_const_iterator_2;
  typedef typename K_VisConstruct_coord_iterator<OldK>        Construct_cartesian_const_iterator_2;
  typedef typename K_VisConstruct_bbox_2<typename OldK::Construct_bbox_2, K_Base>
                                            Construct_bbox_2;


  Construct_point_2
  construct_point_2_object() const
  { return Construct_point_2(); }

  Construct_bbox_2
  construct_bbox_2_object() const
  { return Construct_bbox_2(); }

  Construct_cartesian_const_iterator_2
  construct_cartesian_const_iterator_2_object() const
  { return Construct_cartesian_const_iterator_2(); }

  template < typename Kernel2 >
  struct Base { typedef K_Vis_Base<Kernel2, K_Base>  Type; };
};


template < typename FT_ >
struct K_Vis_Kernel
  : public CGAL::Type_equality_wrapper<
                K_Vis_Base<K_Vis_Kernel<FT_>, CGAL::Cartesian<FT_> >,
                K_Vis_Kernel<FT_> >
{};

#endif // MYKERNEL_H
