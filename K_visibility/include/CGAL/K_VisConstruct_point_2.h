#ifndef K_VISCONSTRUCT_POINT_2_H
#define K_VISCONSTRUCT_POINT_2_H

template <typename K, typename OldK>
class K_VisConstruct_point_2
{
  typedef typename K::RT         RT;
  typedef typename K::Point_2    Point_2;
  typedef typename K::Line_2     Line_2;
  typedef typename Point_2::Rep  Rep;
public:
  typedef Point_2                result_type;
  typedef K_Vis_point_2<OldK> K_Vis_Point ;

  // Note : the CGAL::Return_base_tag is really internal CGAL stuff.
  // Unfortunately it is needed for optimizing away copy-constructions,
  // due to current lack of delegating constructors in the C++ standard.
  Rep // Point_2
  operator()(CGAL::Return_base_tag, CGAL::Origin o) const
  { return Rep(o); }

  Rep // Point_2
  operator()(CGAL::Return_base_tag, const RT& x, const RT& y) const
  { return Rep(x, y); }

  Rep // Point_2
  operator()(CGAL::Return_base_tag, const RT& x, const RT& y, const RT& w) const
  { return Rep(x, y, w); }

  Point_2
  operator()(const CGAL::Origin&) const
  { return K_Vis_Point(0, 0, false); }

  Point_2
  operator()(const RT& x, const RT& y) const
  {
    return K_Vis_Point(x, y, false);
  }

  Point_2
      operator()(const RT& x, const RT& y, bool art) const
  {
      return K_Vis_Point(x, y, art);
  }

  const Point_2&
  operator()(const Point_2 & p) const
  {
    return p;
  }

  Point_2
  operator()(const Line_2& l) const
  {
    typename OldK::Construct_point_2 base_operator;
    Point_2 p = base_operator(l);
    return p;
  }

  Point_2
  operator()(const Line_2& l, int i) const
  {
    typename OldK::Construct_point_2 base_operator;
    return base_operator(l, i);
  }

  // We need this one, as such a functor is in the Filtered_kernel
  Point_2
  operator()(const RT& x, const RT& y, const RT& w) const
  {
    if(w != 1){
      return K_Vis_Point_2(x/w, y/w, 0);
    } else {
      return K_Vis_Point_2(x,y, 0);
    }
  }
};

#endif //K_VISCONSTRUCT_POINT_2_H
