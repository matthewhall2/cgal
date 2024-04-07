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


template<class Kernel>
class K_visibility_region {

    class K_Vis_Point_2 : public Kernel::Point_2 {
    public:

        K_Vis_Point_2(const typename Kernel::FT& x, const typename Kernel::FT& y, bool isArtificial = false) : Kernel::Point_2(x, y), artificial(isArtificial) {
           
        }

        K_Vis_Point_2() : Kernel::Point_2() {

        }

        bool isArtificial() {
            return this.artificial;
        }

        void setPair(K_Vis_Point_2& p1, K_Vis_Point_2& p2) {
            p1.pair = p2;
            p2.pair = p1;
        }

    private:
        bool artificial;
        K_Vis_Point_2 *pair;
    };

    class K_Vis_Kernal : public Kernel {
    public: 
        typedef K_Vis_Point_2 Point_2;
    };


    using Traits = CGAL::Arr_segment_traits_2<K_Vis_Kernal>;
    using Arrangement = CGAL::Arrangement_2<Traits>;
    //typedef Polygon::Vertex_const_iterator EdgeIterator;
    using Point = typename K_Vis_Point_2;
    using Polygon = CGAL::Polygon_2<K_Vis_Kernal>;
    using Transformation = CGAL::Aff_transformation_2<K_Vis_Kernal> ;
    using Segment = CGAL::Segment_2<Traits> ;
    using Line = CGAL::Line_2<K_Vis_Kernal> ;
    using RT = typename Kernel::RT;


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

        std::vector<std::pair<Point *, int>> intersectionPoints;
        int leftIntersectionIndex;
        int rightIntersectionIndex;

        void getUpper();
        void getLower();
        void merge();
        bool isPointHorizontalWithVertex(Point p);
        void projection(Point p);
        void inverseProjection(Point p);
        Polygon lowerProjected;
        Polygon upperProjected;
        void getRadial(Point p);
        boost::unordered_map<std::tuple<RT, RT, RT, RT>, std::vector<Point>> radialIntersectoinList;
        Arrangement upperArr;
};

template<typename Kernel>
K_visibility_region<Kernel>::K_visibility_region(CGAL::Polygon_2<Kernel> p) {
    this->polygon.clear();
    Segment a(Point(0, 0), Point(1, 1));
    Point test(1, 1);
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
}

template<typename Kernel>
void K_visibility_region<Kernel>::getLower() {
    enum  {
        ABOVE=0,
        BELOW=1        
    } state = ABOVE;
    this->lower.p.clear();
    lower.p.push_back(this->polygon.edge(leftIntersectionIndex).start()); // clockwise orientation, start of edge is below line
    std::cout << this->polygon.edge(leftIntersectionIndex).start() << " " << this->polygon.edge(leftIntersectionIndex).end() << "\n";
    lower.p.push_back(*this->intersectionPoints[this->leftIntersectionIndex].first);
    std::cout << *this->intersectionPoints[this->leftIntersectionIndex].first << "\n";
    for (int i = this->leftIntersectionIndex + 1; i < this->intersectionPoints.size(); i++) {
        std::cout << i << "\n";
        if (state == ABOVE) {
            std::cout << "above\n";
            std::cout << this->polygon.edge(i) << "\n";
            if (this->intersectionPoints[i].first == nullptr) {
                continue;
            }
            lower.p.push_back(*this->intersectionPoints[i].first);
            printf("going below\n");
            state = BELOW;
        }
        else if (state == BELOW) {
            std::cout << "below\n";

            lower.p.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                lower.p.push_back(*this->intersectionPoints[i].first);
                state = ABOVE;
            }
        }
    }

    if (leftIntersectionIndex == 0) {
        return;
    }
    
    for (int i = 0; i < this->leftIntersectionIndex; i++) {
        std::cout << i << "\n";
        if (state == ABOVE) {
            std::cout << "above\n";
            if (this->intersectionPoints[i].first == nullptr) {
                continue;
            }
            std::cout << *this->intersectionPoints[i].first << "\n";
            lower.p.push_back(*this->intersectionPoints[i].first);
            state = BELOW;

        }
        else if (state == BELOW) {
            std::cout << "below\n";

            lower.p.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                lower.p.push_back(*this->intersectionPoints[i].first);
                state = ABOVE;
            }
        }
    }
    std::cout << lower.p.edges().size() << "\n";
}

template<typename Kernel>
void K_visibility_region<Kernel>::getUpper() {
    enum {
        ABOVE = 0,
        BELOW = 1
    } state = BELOW;
    this->upper.p.clear();
    upper.p.push_back(this->polygon.edge(rightIntersectionIndex).start()); // clockwise orientation, start of edge is below line
    std::cout << this->polygon.edge(rightIntersectionIndex).start() << " " << this->polygon.edge(rightIntersectionIndex).end() << "\n";
    upper.p.push_back(*this->intersectionPoints[this->rightIntersectionIndex].first);
    std::cout << *this->intersectionPoints[this->rightIntersectionIndex].first << "\n";
    for (int i = this->rightIntersectionIndex + 1; i < this->intersectionPoints.size(); i++) {
        std::cout << i << "\n";
        if (state == BELOW) {
            std::cout << "above\n";
            std::cout << this->polygon.edge(i) << "\n";
            if (this->intersectionPoints[i].first == nullptr) {
                continue;
            }
            upper.p.push_back(*this->intersectionPoints[i].first);
            printf("going below\n");
            state = ABOVE;

        }
        else if (state == ABOVE) {
            std::cout << "below\n";

            upper.p.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                upper.p.push_back(*this->intersectionPoints[i].first);
                state = BELOW;
            }
        }
    }

    if (rightIntersectionIndex == 0) {
        return;
    }

    for (int i = 0; i < this->rightIntersectionIndex; i++) {
        std::cout << i << "\n";
        if (state == BELOW) {
            std::cout << "above\n";
            if (this->intersectionPoints[i].first == nullptr) {
                continue;
            }
            std::cout << *this->intersectionPoints[i].first << "\n";
            upper.p.push_back(*this->intersectionPoints[i].first);
            state = ABOVE;

        }
        else if (state == ABOVE) {
            std::cout << "below\n";

            upper.p.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                upper.p.push_back(*this->intersectionPoints[i].first);
                
                state = BELOW;
            }
        }
    }
   // std::cout << lower.p.edges().size() << "\n";
    //CGAL::draw(upper);
    
}

template<typename Kernel>
typename K_visibility_region<Kernel>::Polygon K_visibility_region<Kernel>::find_visibility_region(int k, Point p) {
    // transformations incase rotation is needed
    this->translateToOrigin = new Transformation(CGAL::TRANSLATION, Vector(-p.x(), -p.y()));
    this->translateBack = new Transformation(CGAL::TRANSLATION, Vector(p.x(), p.y()));
    this->rotate = new Transformation(CGAL::ROTATION, sin(3.1415 / 100), cos(3.1415 / 100));
    if (p.y() == 0) {
        this->moveQueryPoint = new Transformation(CGAL::TRANSLATION, Vector(0, 1));
        p = (*moveQueryPoint)(p);
        this->polygon = (*moveQueryPoint)(this->polygon);
    }

    // rotate polygon if horizontal line intersects a vertex
    while (this->isPointHorizontalWithVertex(p)) {
        for (int i = 0; i < this->polygon.vertices().size(); i++) {
            this->polygon.vertex(i) = (*translateBack)((*rotate)((*translateToOrigin)(polygon.vertex(i))));
        }
    }
    
    
    assert(this->leftIntersectionIndex != -1);
    assert(this->rightIntersectionIndex != -1);

    getUpper();
    getLower();
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
    if ((CGAL::Bounded_side result = this->polygon.bounded_side(p)) == CGAL::ON_BOUNDARY) {
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
    for (int i = 0; i < numEdges; i++) {
        if (p == this->polygon.vertex(i)) {
            return true;
        }
        Segment edge = this->polygon.edge(i);

        const auto intersection = CGAL::intersection(line, edge);
        if (!intersection) {
            this->intersectionPoints.push_back(std::make_pair(nullptr, i));
            continue;
        }

        if (const Segment* s = boost::get<Segment>(&*intersection)) {
            return true;
        }

        const Point* intersectionPoint = boost::get<Point>(&*intersection);
        if (*intersectionPoint == this->polygon.vertex(i)) {
            return true;

        }

        this->intersectionPoints.push_back(std::make_pair(new Point(intersectionPoint->x(), intersectionPoint->y()), i));
        if (*intersectionPoint < left) {
            this->leftIntersectionIndex = this->intersectionPoints.size() - 1;
            left = *intersectionPoint;
        }

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
    this->upperProjected.clear();
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
    CGAL::draw(this->upperProjected);
}

template <class Kernel>
void K_visibility_region<Kernel>::inverseProjection(Point p) {
    //projection is the matrix
    /*1 0 0
    * 0 1 0
    * 0 1 -p.y
    */
    this->upperProjected.clear();
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
    CGAL::draw(this->upperProjected);
}

template <class Kernel>
void K_visibility_region<Kernel>::getRadial(Point p) {
    //this->projection(p);
}
