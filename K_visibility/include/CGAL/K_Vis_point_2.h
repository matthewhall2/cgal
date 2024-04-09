#ifndef K_VIS_POINT_2_H
#define K_VIS_POINT_2_H


#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>

template<typename K_>
class K_Vis_point_2 : public K_::Point_2 {

private:
  
  bool artificial;

public:

  K_Vis_point_2()
    : K_::Point_2(), artificial(false)
  {
    
  }

  
  K_Vis_point_2(const typename K_::FT x, const typename K_::FT y, bool artificial = false) 
    : K_::Point_2(x, y), artificial(artificial)
  {
    
  }

  bool isArtificial() const { return artificial; }

  bool& isArtificial() { return artificial; }


  bool operator==(const K_Vis_point_2 &p) const
  {
    return p.x() == this->x() && p.y() == this->y();
  }

  bool operator!=(const K_Vis_point_2 &p) const
  {
      return !(*this == p);
  }

};

#endif // K_VIS_POINT_2_H
