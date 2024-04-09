#include <CGAL/Polygon_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Kernel_traits.h>
#include <stdbool.h>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <list>
#include <CGAL/Bbox_2.h>
#include <unordered_map>
#include <boost/config.hpp>
#include "K_Vis_Kernel.h"

#include <boost/unordered_map.hpp>
#include <CGAL/Arr_vertical_decomposition_2.h>


#define POINT_TO_DOUBLE(point) std::pair<double, double>((double)point.x(), (double)point.y())
#define EPSILON 10




template<class Kernel>
class K_visibility_region {
    
    using FT = typename Kernel::FT;
    using MK = typename K_Vis_Kernel<FT>;

    using Traits = typename CGAL::Arr_segment_traits_2<MK>;
    using Arrangement = CGAL::Arrangement_2<Traits>;
    //typedef Polygon::Vertex_const_iterator EdgeIterator;
    using Point = typename MK::Point_2;
    using Polygon = CGAL::Polygon_2<MK>;
    using Transformation = CGAL::Aff_transformation_2<MK> ;
    using Segment = typename MK::Segment_2;
    using Line = typename MK::Line_2;
    using Vector = typename CGAL::Vector_2<MK>;

    using Vertex_handle = typename Arrangement::Vertex_handle;
    using Halfedge_handle = typename Arrangement::Halfedge_handle;
    using Face_handle = typename Arrangement::Face_handle;
    using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
    using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
    using Face_const_handle = typename Arrangement::Face_const_handle;

    typedef std::pair<CGAL::Object, CGAL::Object>           Object_pair;
    typedef std::pair<Vertex_const_handle, Object_pair>     Vert_decomp_entry;
    typedef std::list<Vert_decomp_entry>                    Vert_decomp_list;


    public:
        K_visibility_region(CGAL::Polygon_2<Kernel> p);
        Polygon find_visibility_region(int k, CGAL::Point_2<Kernel> p);

    
    private:
        
        Polygon polygon;

        Polygon upper;
        Polygon lower;
        Polygon bounding_box;

        Transformation* rotate;
        Transformation* translateToOrigin;
        Transformation* translateBack;
        Transformation* moveQueryPoint;

        Point queryPoint;
        Point lowestAboveL;
        Point highestBelowL;

        std::vector<std::pair<Point *, int>> intersectionPoints;
        int leftIntersectionIndex;
        int rightIntersectionIndex;

       
        void getLowerUpper();
        void merge();
        bool isPointHorizontalWithVertex(Point p);
        void projection(Point p);
        void inverseProjection(Point p);
        Polygon lowerProjected;
        Polygon upperProjected;
        void getRadial(Point p);
      //  boost::unordered_map<std::tuple<FT, FT, FT, FT>, std::vector<Point>> radialIntersectoinList;
        //std::unordered_map<Point, bool, PointHash> artificialVertexList;
        Arrangement upperArr;
        Arrangement lowerArr;
        std::vector<Segment> upperEdgeList;
        std::vector<Segment> lowerEdgeList;      
};

template<typename Kernel>
K_visibility_region<Kernel>::K_visibility_region(CGAL::Polygon_2<Kernel> p) {
    this->polygon.clear();
    /*Segment a(Point(0, 0), Point(1, 1));
    Segment b(Point(1, 1), Point(2, 1));
    upperEdgeList.push_back(a);
    upperEdgeList.push_back(b);
    insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
    Point test(1, 1);*/
    for (CGAL::Point_2<Kernel> point : p.vertices()) {
     //   Point tt(point.x(), point.y());
        this->polygon.push_back(Point(point.x(), point.y()));
    }

    intersectionPoints.clear();
    if (p.is_empty()) {
        throw std::invalid_argument("Polygon must not be empty.");
    }

    if (!p.is_simple()) {
        throw std::invalid_argument("Polygon must be simple.");
    }

    if (this->polygon.orientation() == CGAL::COUNTERCLOCKWISE) {
        this->polygon.reverse_orientation();
    }
    

    this->lowestAboveL = *this->polygon.top_vertex();
    this->highestBelowL = *this->polygon.bottom_vertex();
    CGAL::draw(p);
}

template<typename Kernel>
void K_visibility_region<Kernel>::getLowerUpper() {
    enum  {
        ABOVE=0,
        BELOW=1        
    } state = ABOVE;

    this->lowerEdgeList.clear();
    FT ly = this->queryPoint.y() - (this->queryPoint.y() - highestBelowL.y()) / EPSILON;
    FT uy = this->queryPoint.y() + (this->highestBelowL.y() - this->queryPoint.y()) / EPSILON;
    Segment e = this->polygon.edge(leftIntersectionIndex);
    Line l = e.supporting_line();
    lowerEdgeList.clear();
    upperEdgeList.clear();
    lowerEdgeList.push_back(Segment(e.start(), Point(l.x_at_y(ly), ly))); // clockwise orientation, start of edge is below line
    upperEdgeList.push_back(Segment(Point(l.x_at_y(uy), uy), e.end()));
   

  //  std::cout << this->polygon.edge(leftIntersectionIndex).start() << " " << this->polygon.edge(leftIntersectionIndex).end() << "\n";
 //   std::cout << *this->intersectionPoints[this->leftIntersectionIndex].first << "\n";
    for (int i = this->leftIntersectionIndex + 1; i < this->intersectionPoints.size(); i++) {
       // printf("point (%d, %d), pointer: %p", this->polygon.vertex(i).x(), this->polygon.vertex(i).y(), this->polygon.vertex(i));
       // std::cout << "uplow: " << this->polygon.vertex(i) <<" "<<  & this->polygon.vertex(i) << "\n";
        e = this->polygon.edge(i);
        l = e.supporting_line();
        std::cout << i << "\n";
        if (state == ABOVE) {
            std::cout << "above\n";
            std::cout << this->polygon.edge(i) << "\n";
            if (this->intersectionPoints[i].first == nullptr) {
               this->upperEdgeList.push_back(e);
                continue;
            }
            Point low(l.x_at_y(ly), ly);
            lowerEdgeList.push_back(Segment(low, e.end()));
           // artificialVertexList[low] = true;

            Point up(l.x_at_y(uy), uy);
            upperEdgeList.push_back(Segment(e.start(), up));
         //   artificialVertexList[up] = true;

            printf("going below\n");
            state = BELOW;
        }
        else if (state == BELOW) {
            std::cout << "below\n";

            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";

                Point low (l.x_at_y(ly), ly);
                lowerEdgeList.push_back(Segment(e.start(), low));
              //  artificialVertexList[low] = true;

                Point up(l.x_at_y(uy), uy);
                upperEdgeList.push_back(Segment(up, e.end()));
              //  artificialVertexList[up] = true;
                state = ABOVE;
            }
            else {
                lowerEdgeList.push_back(e);
            }
        }
    }

    if (leftIntersectionIndex == 0) {
        /*for (const auto& pair : this->artificialVertexList) {
            std::cout << "artificial vertex: (" <<*pair.first <<") " << pair.first << " loc: " << pair.second << "\n";
        }*/
        insert(lowerArr, lowerEdgeList.begin(), lowerEdgeList.end());
        CGAL::draw(lowerArr);
        insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
        CGAL::draw(upperArr);
        Vert_decomp_list vd_list;
       // CGAL::decompose(lowerArr, std::back_inserter(vd_list));
        /*for (auto vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
            const Object_pair& curr = vd_iter->second;
            std::cout << "Vertex (" << vd_iter->first->point() << ") " << &vd_iter->first->point() << " : ";
            Vertex_const_handle vh;
            Halfedge_const_handle hh;
            Face_const_handle fh;
            std::cout << " feature below: ";
            if (CGAL::assign(hh, curr.first)) std::cout << '[' << hh->curve() << ']';
            else if (CGAL::assign(vh, curr.first))
                std::cout << '(' << vh->point() << ')';
            else if (CGAL::assign(fh, curr.first)) std::cout << "NONE";
            else std::cout << "EMPTY";
            std::cout << "   feature above: ";
            if (CGAL::assign(hh, curr.second))
                std::cout << '[' << hh->curve() << "]\n";
            else if (CGAL::assign(vh, curr.second))
                std::cout << '(' << vh->point() << ")\n";
            else if (CGAL::assign(fh, curr.second)) std::cout << "NONE\n";
            else std::cout << "EMPTY\n";
        }*/
        return;
    }
    
    for (int i = 0; i < this->leftIntersectionIndex; i++) {
        e = this->polygon.edge(i);
        l = e.supporting_line();
        std::cout << i << "\n";
        if (state == ABOVE) {
            std::cout << "above\n";
            std::cout << this->polygon.edge(i) << "\n";
            if (this->intersectionPoints[i].first == nullptr) {
                this->upperEdgeList.push_back(e);
                continue;
            }
            Point low(l.x_at_y(ly), ly);
            lowerEdgeList.push_back(Segment(low, e.end()));
            // artificialVertexList[low] = true;

            Point up(l.x_at_y(uy), uy);
            upperEdgeList.push_back(Segment(e.start(), up));
            //   artificialVertexList[up] = true;

            printf("going below\n");
            state = BELOW;
        }
        else if (state == BELOW) {
            std::cout << "below\n";

            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";

                Point low(l.x_at_y(ly), ly);
                lowerEdgeList.push_back(Segment(e.start(), low));
                //  artificialVertexList[low] = true;

                Point up(l.x_at_y(uy), uy);
                upperEdgeList.push_back(Segment(up, e.end()));
                //  artificialVertexList[up] = true;
                state = ABOVE;
            }
            else {
                lowerEdgeList.push_back(e);
            }
        }
    }
    std::cout << "size is " << lowerEdgeList.size() << "\n";
   /* for (const auto& pair : this->artificialVertexList) {
        std::cout << "artificial vertex: (" << *pair.first << ") " << pair.first << " loc: " << pair.second << "\n";
    }*/
    insert(lowerArr, lowerEdgeList.begin(), lowerEdgeList.end());
    CGAL::draw(lowerArr);
    insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
    CGAL::draw(upperArr);
    Vert_decomp_list vd_list;
  //  CGAL::decompose(lowerArr, std::back_inserter(vd_list));
    /*for (auto vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
        const Object_pair& curr = vd_iter->second;
        std::cout << "Vertex (" << vd_iter->first->point() << ") " << &vd_iter->first->point()<<" : ";
        Vertex_const_handle vh;
        Halfedge_const_handle hh;
        Face_const_handle fh;
        std::cout << " feature below: ";
        if (CGAL::assign(hh, curr.first)) std::cout << '[' << hh->curve() << ']';
        else if (CGAL::assign(vh, curr.first))
            std::cout << '(' << vh->point() << ')';
        else if (CGAL::assign(fh, curr.first)) std::cout << "NONE";
        else std::cout << "EMPTY";
        std::cout << "   feature above: ";
        if (CGAL::assign(hh, curr.second))
            std::cout << '[' << hh->curve() << "]\n";
        else if (CGAL::assign(vh, curr.second))
            std::cout << '(' << vh->point() << ")\n";
        else if (CGAL::assign(fh, curr.second)) std::cout << "NONE\n";
        else std::cout << "EMPTY\n";
    }*/
}



template<typename Kernel>
typename K_visibility_region<Kernel>::Polygon K_visibility_region<Kernel>::find_visibility_region(int k, CGAL::Point_2<Kernel> p) {
    // transformations incase rotation is needed
    this->translateToOrigin = new Transformation(CGAL::TRANSLATION, Vector(-p.x(), -p.y()));
    this->translateBack = new Transformation(CGAL::TRANSLATION, Vector(p.x(), p.y()));
    this->rotate = new Transformation(CGAL::ROTATION, sin(3.1415 / 100), cos(3.1415 / 100));
    if (p.y() == 0) {
        this->queryPoint = Point(p.x(), p.y() + 1);
        this->moveQueryPoint = new Transformation(CGAL::TRANSLATION, Vector(0, 1));
      
        this->polygon = CGAL::transform(*this->moveQueryPoint, this->polygon);
    }
    else {
        this->queryPoint = Point(p.x(), p.y());
    }
    auto x = this->polygon.vertex(0);

   //  rotate polygon if horizontal line intersects a vertex
    while (this->isPointHorizontalWithVertex(this->queryPoint)) {
       /* auto test = (*this->translateBack) * (*this->rotate) * (*this->translateToOrigin);
        auto test = CGAL::transform(*this.rotate, this->translateToOrigin);
        this->polygon = CGAL::transform(, this->polygon);
        for (int i = 0; i < this->polygon.vertices().size(); i++) {
            this->polygon.vertex(i) = (*translateBack)((*rotate)((*translateToOrigin)(polygon.vertex(i))));
        }*/
    }
    
    
    assert(this->leftIntersectionIndex != -1);
    assert(this->rightIntersectionIndex != -1);

    getLowerUpper();
    //CGAL::draw(this->lower);
   // CGAL::draw(this->upper);
   // getRadial(p);

    return this->polygon;
}

/*
* checks if horizontal line <l> going through point <p> intersects a vertex of the polygon <P>. 
* If it does not, then list <intersectionPoints> will contain all intersection of line l with the boundary of P in clockwise order
*/
template<typename Kernel>
bool K_visibility_region<Kernel>::isPointHorizontalWithVertex(Point p) {

    // checking if point in inside polygon
    CGAL::Bounded_side result;
    if ((result = this->polygon.bounded_side(p)) == CGAL::ON_BOUNDARY) {
        std::stringstream ss;
        ss << "Point (" << p << ") on boundary. Must be inside polygon";
        throw std::invalid_argument(ss.str());
    }
    else if (result == CGAL::ON_UNBOUNDED_SIDE) {
        std::stringstream ss;
        ss << "Point (" << p << ") outside polygon. Must be inside polygon";
        throw std::invalid_argument(ss.str());
    }

    this->leftIntersectionIndex = -1;
    this->rightIntersectionIndex = -1;
    Point left = p;
    Point right = p;

    Line line(p, Point(p.x() + 1, p.y()));
    int numEdges = this->polygon.edges().size();
    this->intersectionPoints.clear(); // clear before incase of multiple calls due to intersection with vertex
    for (int i = 0; i < numEdges; i++) {
        if (p == this->polygon.vertex(i)) {
            return true;
        }
        Segment edge = this->polygon.edge(i);
        Point vertex = edge.start();

        const auto intersection = CGAL::intersection(line, edge);
        if (!intersection) {
            this->intersectionPoints.push_back(std::make_pair(nullptr, i));
            continue;
        }

        if (vertex.y() < p.y() && vertex.y() > this->highestBelowL.y()) {
            this->highestBelowL = vertex;
        }

        if (vertex.y() > p.y() && vertex.y() < this->lowestAboveL.y()) {
            this->lowestAboveL = vertex;
        }


        // if l is collinear with edge of polygon
        if (const Segment* s = boost::get<Segment>(&*intersection)) {
            return true;
        }

        // if l intersects vertex of polygon
        const Point* intersectionPoint = boost::get<Point>(&*intersection);
        if (*intersectionPoint == this->polygon.vertex(i)) {
            return true;
        }

        // gets leftmost intersection point for use in getLower()
        this->intersectionPoints.push_back(std::make_pair(new Point(intersectionPoint->x(), intersectionPoint->y()), i));
        if (*intersectionPoint < left) {
            this->leftIntersectionIndex = this->intersectionPoints.size() - 1;
            left = *intersectionPoint;
        }

        // gets leftmost intersection point for use in getUpper()
        if (*intersectionPoint > right) {
            this->rightIntersectionIndex = this->intersectionPoints.size() - 1;
            right = *intersectionPoint;
        }
    }

    return false;
}

template <class Kernel>
void K_visibility_region<Kernel>::projection(Point p) {
    //projection is the matrix
    /*1 0 0
    * 0 1 0
    * 0 1 -p.y 
    */
   /* this->upperProjected.clear();
    this->lowerProjected.clear();
    for (Point pp : this->lower.p.vertices()) {
        double x = pp.x();
        double y = pp.y();
        double nz = y + -p.y();
        this->lowerProjected.push_back(Point(x / nz, y / nz));
    }
    CGAL::draw(this->lowerProjected);

    for (Point pp : this->upper.p.p.vertices()) {
        double x = pp.x();
        double y = pp.y();

        double nz = y + -p.y();
        this->upperProjected.push_back(Point(x / nz, y / nz));
    }
    CGAL::draw(this->upperProjected);*/
}

template <class Kernel>
void K_visibility_region<Kernel>::inverseProjection(Point p) {
    //projection is the matrix
    /*1 0 0
    * 0 1 0
    * 0 1 -p.y
    */
    /*this->upperProjected.clear();
    this->lowerProjected.clear();
    for (Point pp : this->lower.p.vertices()) {
        double x = pp.x();
        double y = pp.y();
        double nz = y + -p.y();
        this->lowerProjected.push_back(Point(x / nz, y / nz));
    }
    CGAL::draw(this->lowerProjected);

    for (Point pp : this->upper.p.p.vertices()) {
        double x = pp.x();
        double y = pp.y();

        double nz = y + -p.y();
        this->upperProjected.push_back(Point(x / nz, y / nz));
    }
    CGAL::draw(this->upperProjected);*/
}

template <class Kernel>
void K_visibility_region<Kernel>::getRadial(Point p) {
    //this->projection(p);
}
