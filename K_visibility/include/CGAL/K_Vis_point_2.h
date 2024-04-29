#ifndef K_VIS_POINT_2_H
#define K_VIS_POINT_2_H


#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>

template<typename K_>
class K_Vis_point_2 : public K_::Point_2 {
private:
  
  bool artificial;
  int id_;
  bool isTip_;
  typename K_::FT testx;
  typename K_::FT testy;
public:

  K_Vis_point_2(): K_::Point_2(), artificial(false), id_(-1), isTip_(false)
  {
     
  }

  
  K_Vis_point_2(const typename K_::FT x, const typename K_::FT y, bool artificial = false, int id=-1) 
    : K_::Point_2(x, y), testx(x), testy(y), artificial(artificial), id_(id), isTip_(false)
  {
     
  }

  K_Vis_point_2(const typename K_::FT x, const typename K_::FT y, const typename K_::FT w,  bool artificial = false, int id = -1)
      : K_::Point_2(x, y, w), testx(x), testy(y), artificial(artificial), id_(id)
  {

  }

  std::string toString() {
      std::string s = "";
      std::ostringstream  o;
   o <<  "[(" << this->hx() << " " << this->hy() << " " << this->hw() << ")" << " id: " << this->id() << " artificial: " << (this->isArtificial() ? "true" : "false") << "]";
   return o.str();
  }

  bool isArtificial() const { return artificial; }

  bool& isArtificial() { return artificial; }

  int id() const { return id_; }

  int& id() { return id_; }

  bool isTip() const { return isTip_; }
  bool& isTip() { return isTip_; }

  void setId(int id) {
      id_ = id;
  }


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
