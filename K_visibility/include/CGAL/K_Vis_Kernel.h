#ifndef MYKERNEL_H
#define MYKERNEL_H

#include <CGAL/Cartesian.h>
#include "K_Vis_point_2.h"
#include "K_VisSegment_2.h"
#include "K_VisConstruct_bbox_2.h"
#include "K_VisConstruct_coord_iterator.h"
#include "K_VisConstruct_point_2.h"
//#include "K_VisSeg_2.h"
#include "K_VisConstructSegment_2.h"
//#include "K_VisConstructSeg_2.h"

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
  typedef typename const K_Base::FT* Cartesian_const_iterator_2;
  typedef typename K_VisConstruct_coord_iterator<OldK>        Construct_cartesian_const_iterator_2;
 typedef typename K_VisConstruct_bbox_2<typename OldK::Construct_bbox_2, K_Base>
                                            Construct_bbox_2;
 typedef typename K_VisConstructSegment_2<Kernel> Construct_segment_2;


  Construct_point_2
  construct_point_2_object() const
  {
 return Construct_point_2();
}

Construct_bbox_2
construct_bbox_2_object() const
{
return Construct_bbox_2();
}

Construct_cartesian_const_iterator_2
construct_cartesian_const_iterator_2_object() const
{
return Construct_cartesian_const_iterator_2();
}

Construct_segment_2
    construct_segment_2_object() const
{
    return Construct_segment_2();
}

template < typename Kernel2 >
struct Base { typedef K_Vis_Base<Kernel2, K_Base>  Type; };
};


//template < typename K_, typename K_Base >
//class K_Vis_Base2
//    : public K_Base::template Base<K_>::Type
//{
//  typedef typename K_Base::template Base<K_>::Type   OldK;
//public:
//
//  typedef typename K_                                Kernel;
//  typedef typename K_VisSeg_2<OldK>               Segment_2;
//
//  typedef typename K_VisConstructSeg_2<Kernel, OldK> Construct_segment_2;
//
//
//Construct_segment_2
//    construct_segment_2_object() const
//{
//    return Construct_segment_2();
//}
//
//template < typename Kernel2 >
//struct Base { typedef K_Vis_Base2<Kernel2, K_Base>  Type; };
//};
//
template < typename FT_ >
struct K_Vis_Kernel
    : public CGAL::Type_equality_wrapper<
    K_Vis_Base<K_Vis_Kernel<FT_>, CGAL::Cartesian<FT_> >,
    K_Vis_Kernel<FT_> >
{};
//
//namespace CGAL {
//    template < class R >
//    std::ostringstream&
//        operator<<(std::ostringstream& os, const CGAL::Point_2<R>& s)
//    {
//        switch (CGAL::IO::get_mode(os)) {
//        case CGAL::IO::ASCII:
//            return os "[(" << s.hx() << " " << s.hy() << " " << s.hw() << ")" << " id: " << s.id() << " artificial: " << (s.isArtificial() ? "true" : "false") << "]";
//        case CGAL::IO::BINARY:
//            return os << s.source() << s.target() << " id: " << s.id();
//        default:
//            return os << "K_Vis_point_2(" << s.source() << ", " << s.target() << ")";
//        }
//    }
//}
//
//template < typename FT_ >
//struct K_Vis_Kernel2
//    : public CGAL::Type_equality_wrapper<
//    K_Vis_Base2<K_Vis_Kernel2<FT_>, K_Vis_Kernel<FT_> >,
//    K_Vis_Kernel2<FT_> >
//{};

#endif // MYKERNEL_H

template <typename K>
class K_Vis_Transformation_2 : public CGAL::Aff_transformation_2<K> {
public:
    using Point_2 = typename K::Point_2;
    using CGAL::Aff_transformation_2<K>::Aff_transformation_2;
    Point_2
        operator()(const Point_2& p) const
    {
        Point_2 r = transform(p);
        r.id() = p.id();
        r.isArtificial() = p.isArtificial();
        r.isTip() = p.isTip();
        r.isReflex() = p.isReflex();
        r.baseAwayFromQP() = p.baseAwayFromQP();
        return r;
    }
};
