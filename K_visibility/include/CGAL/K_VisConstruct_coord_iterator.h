#ifndef K_VisCONSTRUCT_COORD_ITERATOR_H
#define K_VisCONSTRUCT_COORD_ITERATOR_H

template <typename K_>
class K_VisConstruct_coord_iterator {
public:
  const typename K_::FT* operator()(const K_Vis_point_2<K_>& p)
  {
    return &p.x();
  }

  const typename K_::FT* operator()(const K_Vis_point_2<K_>& p, int)
  {
    const K_::FT* pyptr = &p.y();
    pyptr++;
    return pyptr;
  }
};

#endif //K_VisCONSTRUCT_COORD_ITERATOR_H
