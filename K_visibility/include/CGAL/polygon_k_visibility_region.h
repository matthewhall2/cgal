#include <CGAL/Polygon_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/intersections.h>
#include <stdbool.h>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <list>



template<class Kernel>
class K_visibility_region {
    typedef CGAL::Polygon_2<Kernel> Polygon;
    //typedef Polygon::Vertex_const_iterator EdgeIterator;
    typedef CGAL::Point_2<Kernel> Point;
    typedef CGAL::Aff_transformation_2<Kernel> Transformation;
    typedef CGAL::Segment_2<Kernel> Segment;
    typedef CGAL::Line_2<Kernel> Line;


    public:
        K_visibility_region(Polygon p);
        Polygon find_visibility_region(int k, Point p);

    
    private:
        Polygon polygon;

        Polygon upper;
        Polygon lower;
        Polygon bounding_box;

        Transformation* rotate;
        Transformation* translateToOrigin;
        Transformation* translateBack;

        std::vector<std::pair<Point *, int>> intersectionPoints;
        int leftIntersectionIndex;
        int rightIntersectionIndex;

        void getUpper();
        void getLower();
        void merge();
        bool isPointHorizontalWithVertex(Point p);
};

template<typename Kernel>
K_visibility_region<Kernel>::K_visibility_region(Polygon p) : polygon(p) {
    intersectionPoints.clear();
    if (this->polygon.orientation() == CGAL::COUNTERCLOCKWISE) {
        this->polygon.reverse_orientation();
    }

    assert(this->polygon.orientation() == CGAL::CLOCKWISE);
}

template<typename Kernel>
void K_visibility_region<Kernel>::getLower() {
    enum  {
        ABOVE=0,
        BELOW=1        
    } state = ABOVE;
    this->lower.clear();
    lower.push_back(this->polygon.edge(leftIntersectionIndex).start()); // clockwise orientation, start of edge is below line
    std::cout << this->polygon.edge(leftIntersectionIndex).start() << " " << this->polygon.edge(leftIntersectionIndex).end() << "\n";
    lower.push_back(*this->intersectionPoints[this->leftIntersectionIndex].first);
    std::cout << *this->intersectionPoints[this->leftIntersectionIndex].first << "\n";
    for (int i = this->leftIntersectionIndex + 1; i < this->intersectionPoints.size(); i++) {
        std::cout << i << "\n";
        if (state == ABOVE) {
            std::cout << "above\n";
            std::cout << this->polygon.edge(i) << "\n";
            if (this->intersectionPoints[i].first == nullptr) {
                continue;
            }
            lower.push_back(*this->intersectionPoints[i].first);
            printf("going below\n");
            state = BELOW;

        }
        else if (state == BELOW) {
            std::cout << "below\n";

            lower.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                lower.push_back(*this->intersectionPoints[i].first);
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
            lower.push_back(*this->intersectionPoints[i].first);
            state = BELOW;

        }
        else if (state == BELOW) {
            std::cout << "below\n";

            lower.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                lower.push_back(*this->intersectionPoints[i].first);
                state = ABOVE;
            }
        }
    }
    std::cout << lower.edges().size() << "\n";
}

template<typename Kernel>
void K_visibility_region<Kernel>::getUpper() {
    enum {
        ABOVE = 0,
        BELOW = 1
    } state = BELOW;
    this->upper.clear();
    upper.push_back(this->polygon.edge(rightIntersectionIndex).start()); // clockwise orientation, start of edge is below line
    std::cout << this->polygon.edge(rightIntersectionIndex).start() << " " << this->polygon.edge(rightIntersectionIndex).end() << "\n";
    upper.push_back(*this->intersectionPoints[this->rightIntersectionIndex].first);
    std::cout << *this->intersectionPoints[this->rightIntersectionIndex].first << "\n";
    for (int i = this->rightIntersectionIndex + 1; i < this->intersectionPoints.size(); i++) {
        std::cout << i << "\n";
        if (state == BELOW) {
            std::cout << "above\n";
            std::cout << this->polygon.edge(i) << "\n";
            if (this->intersectionPoints[i].first == nullptr) {
                continue;
            }
            upper.push_back(*this->intersectionPoints[i].first);
            printf("going below\n");
            state = ABOVE;

        }
        else if (state == ABOVE) {
            std::cout << "below\n";

            upper.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                upper.push_back(*this->intersectionPoints[i].first);
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
            upper.push_back(*this->intersectionPoints[i].first);
            state = ABOVE;

        }
        else if (state == ABOVE) {
            std::cout << "below\n";

            upper.push_back(this->polygon.vertex(i));
            if (this->intersectionPoints[i].first != nullptr) {
                std::cout << *this->intersectionPoints[i].first << "\n";
                upper.push_back(*this->intersectionPoints[i].first);
                state = BELOW;
            }
        }
    }
    std::cout << lower.edges().size() << "\n";
    //CGAL::draw(upper);
    
}

template<typename Kernel>
CGAL::Polygon_2<Kernel> K_visibility_region<Kernel>::find_visibility_region(int k, CGAL::Point_2<Kernel> p) {
    // transformations incase rotation is needed
    this->translateToOrigin = new Transformation(CGAL::TRANSLATION, Vector(-p.x(), -p.y()));
    this->translateBack = new Transformation(CGAL::TRANSLATION, Vector(p.x(), p.y()));
    this->rotate = new Transformation(CGAL::ROTATION, sin(3.1415 / 100), cos(3.1415 / 100));

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
    CGAL::draw(this->lower);
    CGAL::draw(this->upper);



    return this->polygon;
}

template<typename Kernel>
bool K_visibility_region<Kernel>::isPointHorizontalWithVertex(CGAL::Point_2<Kernel> p) {
    if (p.x() < polygon.left_vertex()->x() || p.x() > polygon.right_vertex()->x() || p.y() < polygon.bottom_vertex()->y() || p.y() > polygon.top_vertex()->y()) {
        std::stringstream ss;
        ss << "Point (" << p << ") not in polygon";
        throw std::invalid_argument(ss.str());
    }

    this->leftIntersectionIndex = -1;
    this->rightIntersectionIndex = -1;
    Point left = p;
    Point right = p;

    Line line(Segment(p, Point_2(p.hx() + 1, p.hy())));
    int numEdges = this->polygon.edges().size();
    for (int i = 0; i < numEdges; i++) {
        Segment edge = this->polygon.edge(i);

        const auto intersection = CGAL::intersection(line, edge);
        if (!intersection) {
            this->intersectionPoints.push_back(std::make_pair(nullptr, i));
            continue;
        }

        if (const Segment_2* s = boost::get<Segment_2>(&*intersection)) {
            return true;
        }

        const Point* intersectionPoint = boost::get<Point>(&*intersection);
        if (*intersectionPoint == this->polygon.vertex(i)) {
            return true;

        }

        this->intersectionPoints.push_back(std::make_pair(new Point_2(intersectionPoint->x(), intersectionPoint->y()), i));
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
