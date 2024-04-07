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
#include <boost/unordered_map.hpp>
#include <CGAL/Arr_vertical_decomposition_2.h>


#define EPSILON 10




template<class Kernel>
class K_visibility_region {

    


    using Traits = CGAL::Arr_segment_traits_2<Kernel>;
    using Arrangement = CGAL::Arrangement_2<Traits>;
    //typedef Polygon::Vertex_const_iterator EdgeIterator;
    using Point = CGAL::Point_2<Kernel>;
    using Polygon = CGAL::Polygon_2<Kernel>;
    using Transformation = CGAL::Aff_transformation_2<Kernel> ;
    using Segment = CGAL::Segment_2<Kernel> ;
    using Line = CGAL::Line_2<Kernel> ;
    using RT = typename Kernel::RT;

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
        Polygon find_visibility_region(int k, Point p);

    
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
        boost::unordered_map<std::tuple<RT, RT, RT, RT>, std::vector<Point>> radialIntersectoinList;
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
    this->lowestAboveL = *p.top_vertex();
    this->highestBelowL = *p.bottom_vertex();
    CGAL::draw(p);
}

template<typename Kernel>
void K_visibility_region<Kernel>::getLowerUpper() {
    enum  {
        ABOVE=0,
        BELOW=1        
    } state = ABOVE;

    this->lowerEdgeList.clear();
    RT ly = this->queryPoint.y() - (this->queryPoint.y() - highestBelowL.y()) / EPSILON;
    RT uy = this->queryPoint.y() + (this->highestBelowL.y() - this->queryPoint.y()) / EPSILON;
    Segment e = this->polygon.edge(leftIntersectionIndex);
    Line l = e.supporting_line();
    lowerEdgeList.clear();
    upperEdgeList.clear();
    lowerEdgeList.push_back(Segment(e.start(), Point(l.x_at_y(ly), ly))); // clockwise orientation, start of edge is below line
    upperEdgeList.push_back(Segment(Point(l.x_at_y(uy), uy), e.end()));
   

    std::cout << this->polygon.edge(leftIntersectionIndex).start() << " " << this->polygon.edge(leftIntersectionIndex).end() << "\n";
    std::cout << *this->intersectionPoints[this->leftIntersectionIndex].first << "\n";
    for (int i = this->leftIntersectionIndex + 1; i < this->intersectionPoints.size(); i++) {
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
            lowerEdgeList.push_back(Segment(Point(l.x_at_y(ly), ly), e.end()));
            upperEdgeList.push_back(Segment(e.start(), Point(l.x_at_y(uy), uy)));
            printf("going below\n");
            state = BELOW;
        }
        else if (state == BELOW) {
            std::cout << "below\n";

            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                lowerEdgeList.push_back(Segment(e.start(), Point(l.x_at_y(ly), ly)));
                upperEdgeList.push_back(Segment(Point(l.x_at_y(uy), uy), e.end()));
                state = ABOVE;
            }
            else {
                lowerEdgeList.push_back(e);
            }
        }
    }

    if (leftIntersectionIndex == 0) {
        insert(lowerArr, lowerEdgeList.begin(), lowerEdgeList.end());
        CGAL::draw(lowerArr);
        insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
        CGAL::draw(upperArr);
        Vert_decomp_list vd_list;
        CGAL::decompose(lowerArr, std::back_inserter(vd_list));

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
            lowerEdgeList.push_back(Segment(Point(l.x_at_y(ly), ly), e.end()));
            upperEdgeList.push_back(Segment(e.start(), Point(l.x_at_y(uy), uy)));
            printf("going below\n");
            state = BELOW;
        }
        else if (state == BELOW) {
            std::cout << "below\n";

            // lowerEdgeList.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                lowerEdgeList.push_back(Segment(e.start(), Point(l.x_at_y(ly), ly)));
                upperEdgeList.push_back(Segment(Point(l.x_at_y(uy), uy), e.end()));

                state = ABOVE;
            }
            else {
                lowerEdgeList.push_back(e);
            }
        }
    }
    std::cout << "size is " << lowerEdgeList.size() << "\n";
    insert(lowerArr, lowerEdgeList.begin(), lowerEdgeList.end());
    CGAL::draw(lowerArr);
    insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
    CGAL::draw(upperArr);
}



template<typename Kernel>
typename K_visibility_region<Kernel>::Polygon K_visibility_region<Kernel>::find_visibility_region(int k, Point p) {
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
        this->queryPoint = p;
    }

    // rotate polygon if horizontal line intersects a vertex
    while (this->isPointHorizontalWithVertex(p)) {
        for (int i = 0; i < this->polygon.vertices().size(); i++) {
            this->polygon.vertex(i) = (*translateBack)((*rotate)((*translateToOrigin)(polygon.vertex(i))));
        }
    }
    
    
    assert(this->leftIntersectionIndex != -1);
    assert(this->rightIntersectionIndex != -1);

    getLowerUpper();
    //CGAL::draw(this->lower);
   // CGAL::draw(this->upper);
    getRadial(p);

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
