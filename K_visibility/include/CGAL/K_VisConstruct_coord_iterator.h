#ifndef K_VisCONSTRUCT_COORD_ITERATOR_H
#define K_VisCONSTRUCT_COORD_ITERATOR_H

template <typename K_>
class K_VisConstruct_coord_iterator {
public:
  const double* operator()(const K_Vis_point_2<K_>& p)
  {
    return &p.x();
  }

  const double* operator()(const K_Vis_point_2<K_>& p, int)
  {
    const double* pyptr = &p.y();
    pyptr++;
    return pyptr;
  }
};

#endif //K_VisCONSTRUCT_COORD_ITERATOR_H
