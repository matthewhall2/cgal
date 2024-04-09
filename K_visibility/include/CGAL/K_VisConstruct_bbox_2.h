#ifndef K_VISCONSTRUCT_BBOX_2_H
#define K_VISCONSTRUCT_BBOX_2_H


template <class ConstructBbox_2, typename K_>
class K_VisConstruct_bbox_2 : public ConstructBbox_2 {
public:
  using ConstructBbox_2::operator();

  CGAL::Bbox_2 operator()(const K_Vis_point_2<K_>& p) const {
    return CGAL::Bbox_2(p.x(), p.y(), p.x(), p.y());
  }
};

#endif //K_VISCONSTRUCT_BBOX_2_H
