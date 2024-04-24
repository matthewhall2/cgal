#include <CGAL/draw_polygon_set_2.h>
#include <CGAL/Aff_transformation_2.h>
template<typename K_>
class K_VisSeg_2 : public K_VisSegment_2<K_> {
    private:
        int id_;

    public:
        

        

        int& id() { return id_; }

        int id() const { return id_; }


};
