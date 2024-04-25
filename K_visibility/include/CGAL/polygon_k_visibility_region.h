#include <CGAL/Polygon_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/draw_polygon_set_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/draw_arrangement_2.h>
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
#include "K_VisArr_segment_traits_2.h"
#include <CGAL/Unique_hash_map.h>
#include <boost/unordered_map.hpp>
#include <CGAL/Arr_vertical_decomposition_2.h>


#define EPSILON 1/100



    template<class Kernel>
    class K_visibility_region {

        using FT = typename Kernel::FT;
        using RT = typename Kernel::RT;
        using MK = typename K_Vis_Kernel<FT>;

        using Traits = typename CGAL::Arr_segment_traits_2<MK>;
        using GPS_Traits = typename CGAL::Gps_segment_traits_2<MK>;
        using Arrangement = CGAL::Arrangement_2<Traits>;
        //  using Poly_traits = typename CGAL::Arr_polyline_traits_2<MK>;
        //  using Poly_arr = CGAL::Arrangement_2<Poly_traits>;
          //typedef Polygon::Vertex_const_iterator EdgeIterator;
        using Point = typename MK::Point_2;
        using Polygon = CGAL::Polygon_2<MK>;
        using Transformation = K_Vis_Transformation_2<MK>;
        using Segment = typename MK::Segment_2;
        using Line = typename Traits::Line_2;
        using Vector = typename CGAL::Vector_2<MK>;
        using Polygon_set = typename CGAL::Polygon_set_2<MK>;

        using Vertex_handle = typename Arrangement::Vertex_handle;
        using Halfedge_handle = typename Arrangement::Halfedge_handle;
        using Face_handle = typename Arrangement::Face_handle;
        using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
        using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
        using Face_const_handle = typename Arrangement::Face_const_handle;
        using Face_iterator = typename Arrangement::Face_iterator;
        using Face_handle = typename Arrangement::Face_handle;
        using Halfedge_iterator = typename Arrangement::Halfedge_iterator;
        using Halfedge_handle = typename Arrangement::Halfedge_handle;



        typedef std::pair<CGAL::Object, CGAL::Object>           Object_pair;
        typedef std::pair<Vertex_const_handle, Object_pair>     Vert_decomp_entry;
        typedef std::list<Vert_decomp_entry>                    Vert_decomp_list;


    public:
        K_visibility_region(CGAL::Polygon_2<Kernel> p);
        Polygon find_visibility_region(int k, CGAL::Point_2<Kernel> p);


    private:

        Polygon polygon;

        Polygon bounding_box;

        Transformation* rotate;
        Transformation* translateToOrigin;
        Transformation* translateBack;
        Transformation* moveQueryPoint;

        Point queryPoint;
        Point lowestAboveL;
        Point highestBelowL;

        std::vector<std::pair<Point*, int>> intersectionPoints;
        int leftIntersectionIndex;
        int rightIntersectionIndex;
        Face_handle unboundedFace;

        void getLowerUpper();
        void lowerUpperHelper(int start, int end, int&, FT ly, FT uy, std::vector<Segment>&, std::vector<Segment>&);
        void merge();
        bool isPointHorizontalWithVertex(Point p);
        Point projection(Point p) {
            // projection is the matrix
            /*1 0 0
            * 0 1 0
            * 0 1 -qp.y
            */
            //  std::cout << "projecting point: " << p << "\n";
            std::cout << "projecting with id: " << p.id() << std::endl;
            FT x = p.x();
            FT y = p.y();
            FT nz = y - this->queryPoint.y();
            //  std::cout << "query point: " << this->queryPoint << " hw = " << nz <<  std::endl;
            this->hws[p.id()] = nz;
            //    std::cout << "proj dif: " << nz << "\n";
            Point P(x / nz, y / nz, 1);
            //auto P = new Point(x , y , nz);
            P.isArtificial() = p.isArtificial();
            P.id() = p.id();
           // assert(hws.find(-1) == hws.end());
            return  P;
        }

        Point inverseProjection(Point p) {
            // projection is the matrix
           /*1 0 0
           * 0 1 0
           * 0 1/qp.y -1/qp.y
           */
            std::cout << "inv proj id" << p.id() << std::endl;
            std::cout << "inv projecting point: " << p << "\n";
            std::cout << "homo: " << p.hx() << " " << p.hy() << std::endl;
            FT x = p.x();
            FT y = p.y();
            //  FT nz = y - this->queryPoint.y();
            FT nz = CGAL::inverse(this->queryPoint.y()) * y - CGAL::inverse(this->queryPoint.y()) * this->hws[p.id()];
            std::cout << "inv proj query point: " << this->queryPoint << " hw = " << this->hws[p.id()] << std::endl;

            //     std::cout << "proj dif: " << nz << "\n";
                 //std::cout << "query point: " << this->queryPoint << std::endl;
            Point P(x * this->hws[p.id()], y * this->hws[p.id()], 1);
            //auto P = new Point(p.x() * p.hw(), p.y() * p.hw(), 1);
            P.isArtificial() = p.isArtificial();
            return  P;
        }

        Polygon lowerProjected;
        Polygon upperProjected;
        void getRadial(Point p);
        void radialHelper(Point p, Arrangement arr, int &artCounter);
        //    CGAL::Unique_hash_map<std::tuple<FT, FT, FT, FT>, bool> ttt;
        //   boost::unordered_map<std::tuple<FT, FT, FT, FT>, std::vector<Point>> radialIntersectoinList;
          //  std::unordered_map<Point, bool> artificialVertexList;
        Arrangement upperArr;
        Arrangement lowerArr;
        Arrangement arr;
        std::vector<Segment> upperEdgeList;
        std::vector<Segment> lowerEdgeList;
        std::vector<Segment> radialList;
        std::unordered_map<int, FT> hws;

        void projection_insert(std::vector<Segment>& vec, const Segment& seg);

        
        Point interpolate(Halfedge_const_handle &hh, Vertex_const_handle &vh) {
            Point source = hh->source()->point();
            Point target = hh->target()->point();
            Point p = vh->point();
            FT lineLength = CGAL::squared_distance(source, target);
            FT distFromSource = CGAL::squared_distance(source, p);
            FT t = distFromSource / lineLength;
            assert(this->hws.find(source.id()) != hws.end());
            assert(this->hws.find(target.id()) != hws.end());

            FT sourceW = this->hws[source.id()];
            assert(sourceW != 0);
            FT targetW = this->hws[target.id()];
            assert(targetW != 0);
            FT w = sourceW + t * (targetW - sourceW);
            Line l = (Line)hh->curve();
            return Point(p.x() * w, l.y_at_x(p.x()) * w, 1);
        }


        void inv_projection_insert(const Segment& seg);
        void inv_projection_insert(Halfedge_const_handle& bhh, Halfedge_const_handle& uhh, Vertex_const_handle &qvh);
        void inv_projection_insert(Vertex_const_handle &bvh, Vertex_const_handle &uvh, Vertex_const_handle &);

        void inv_projection_insert(Halfedge_const_handle &bhh, Vertex_const_handle &uhh, Vertex_const_handle &);
        void inv_projection_insert(Vertex_const_handle &bhh, Halfedge_const_handle &uhh, Vertex_const_handle &);

        void inv_projection_insert(Vertex_const_handle& bhh,Vertex_const_handle&);
        void inv_projection_insert(Halfedge_const_handle& uhh, Vertex_const_handle&);



    };

    template<typename Kernel>
    K_visibility_region<Kernel>::K_visibility_region(CGAL::Polygon_2<Kernel> p) {
        this->polygon.clear();
        /*Segment a(Point(0, 0), Point(1, 1));
        Segment b(Point(1, 1), Point(2, 1));
        projection_insert(upperEdgeList, a);
        projection_insert(upperEdgeList, b);
        insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
        Point test(1, 1);*/
        int i = 0;
        for (CGAL::Point_2<Kernel> point : p.vertices()) {
            //   Point tt(point.x(), point.y());
            Point p(point.x(), point.y());
            p.id() = i;
            this->polygon.push_back(p);
            i += 1;
        }

        for (Point point : this->polygon.vertices()) {
            //   Point tt(point.x(), point.y());
            std::cout << point.id() << std::endl;
        }

        int testSegId = 1000;
        for (auto edge = this->polygon.edges_begin(); edge != this->polygon.edges_end(); edge++) {            //   Point tt(point.x(), point.y());
            std::cout << edge->toString() << std::endl;
            edge->id() = testSegId++;
        }


        std::vector<Segment> x;
        x.push_back(Segment(Point(0, 0), Point(1, 1)));
        x.push_back(Segment(Point(2, 2), Point(2, 1)));
        Arrangement testa;
        insert(testa, x.begin(), x.end());
        // CGAL::draw(testa);

        Vert_decomp_list vd_list;
        CGAL::decompose(testa, std::back_inserter(vd_list));
        int count = 0;
        for (auto vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
            count += 1;
        }
        std::cout << "count in test: " << count << std::endl;


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
        //  CGAL::draw(p);
    }

    template<typename Kernel>
    void K_visibility_region<Kernel>::lowerUpperHelper(int start, int end, int &artCounter, FT ly, FT uy, std::vector<Segment> &testLower, std::vector<Segment> &testUpper) {
        static enum {
            ABOVE = 0,
            BELOW = 1
        } state = ABOVE;
        Segment e;
        Line l;
        int nv = this->polygon.vertices().size();

        for (int i = start; i < end; i++) {
            e = *(this->polygon.edges_begin() + i);
            std::cout << "edge has id "  << e.id() << std::endl;
          //  Point p1 = (this->polygon[i]);
          //  Point p2 = (this->polygon[i + 1]);

            std::cout << "dividing: edge has ids " << e[0].id() << " " << e.target().id() << " artificial?: " << (e.source().isArtificial() ? "true" : "false") << " " << (e.target().isArtificial() ? "true" : "false") << std::endl;
          //  std::cout << "by point: edge has ids" << p1.id() << " " << p2.id() << " artificial?: " << (p1.isArtificial() ? "true" : "false") << " " << (p2.isArtificial() ? "true" : "false") << std::endl;

            l = e.supporting_line();
            if (state == ABOVE) {
                std::cout << this->polygon.edge(i) << "\n";
                if (this->intersectionPoints[i].first == nullptr) {
                    projection_insert(upperEdgeList, e);
                    testUpper.push_back(e);
                    continue;
                }
                Point low(l.x_at_y(ly), ly);
                low.id() = nv + artCounter++;
                low.isArtificial() = true;
                Segment lowSeg(low, e.end());
                // its only inserting artificial vertices in the list? only articifial have valid ids?
                assert((lowSeg.target().id() == nv + artCounter - 1) || (lowSeg.source().id() == nv + artCounter - 1));

                lowSeg.id() = e.id();
                projection_insert(lowerEdgeList, lowSeg);
                testLower.push_back(lowSeg);

                Point up(l.x_at_y(uy), uy);
                up.isArtificial() = true;
                up.id() = nv + artCounter++;
                Segment upSeg(e.start(), up);
                assert((upSeg.target().id() == nv + artCounter - 1) || (upSeg.source().id() == nv + artCounter - 1));

                upSeg.id() = e.id();
                projection_insert(upperEdgeList, upSeg);
                testUpper.push_back(upSeg);

                state = BELOW;
            }
            else if (state == BELOW) {

                if (this->intersectionPoints[i].first != nullptr) {
                    std::cout << *this->intersectionPoints[i].first << "\n";

                    Point low(l.x_at_y(ly), ly);
                    low.isArtificial() = true;
                    low.id() = nv + artCounter++;
                    Segment lowSeg(e.start(), low);
                    assert((lowSeg.target().id() == nv + artCounter - 1) || (lowSeg.source().id() == nv + artCounter - 1));
                    lowSeg.id() = e.id();
                    projection_insert(lowerEdgeList, lowSeg);
                    testLower.push_back(lowSeg);

                    Point up(l.x_at_y(uy), uy);
                    up.isArtificial() = true;
                    up.id() = nv + artCounter++;

                    Segment upSeg(up, e.end());
                    upSeg.id() = e.id();
                    assert((upSeg.target().id() == nv + artCounter - 1) || (upSeg.source().id() == nv + artCounter - 1));

                    projection_insert(upperEdgeList, upSeg);
                    testUpper.push_back(upSeg);
                    state = ABOVE;
                }
                else {
                    projection_insert(lowerEdgeList, e);
                    testLower.push_back(e);
                }
            }
        }
    }

    template<typename Kernel>
    void K_visibility_region<Kernel>::getLowerUpper() {
        
        std::vector<Segment> testLower;
        std::vector<Segment> testUpper;
        this->lowerEdgeList.clear();
        int artCounter = 0;
        int nv = this->polygon.vertices().size();
        assert(this->queryPoint.y() - highestBelowL.y() > 0);
        assert(lowestAboveL.y() - this->queryPoint.y() > 0);

        FT ly = this->queryPoint.y() - (this->queryPoint.y() - highestBelowL.y()) * EPSILON;

        FT uy = this->queryPoint.y() + (this->lowestAboveL.y() - this->queryPoint.y()) * EPSILON;
        Segment e = this->polygon.edge(leftIntersectionIndex);
        Line l = e.supporting_line();
        lowerEdgeList.clear();
        upperEdgeList.clear();
        Point low(l.x_at_y(ly), ly);
        low.id() = nv + artCounter++;
        low.isArtificial() = true;
        Segment lowSeg(e.start(), low);
        projection_insert(lowerEdgeList, lowSeg); // clockwise orientation, start of edge is below line

        Point up(l.x_at_y(uy), uy);
        up.isArtificial() = true;
        up.id() = nv + artCounter++;
        Segment upSeg(up, e.end());
        projection_insert(upperEdgeList, upSeg);
        //   ttt[std::make_tuple(lowSeg.start().x(), lowSeg.start().y(), upSeg.end().x(), upSeg.end().y())] = true;
        this->lowerUpperHelper(this->leftIntersectionIndex + 1, this->intersectionPoints.size(), artCounter, ly, uy, testLower, testUpper);
        this->lowerUpperHelper(0, this->leftIntersectionIndex + 1, artCounter, ly, uy, testLower, testUpper);

        std::cout << "trgethsrh";
        std::cout << "size is " << lowerEdgeList.size() << std::endl;
        std::cout << "number of homogeneous: " << this->hws.size() << std::endl;
        Arrangement testArrLower;
        insert(testArrLower, testLower.begin(), testLower.end());
        std::cout << "There are " << artCounter << " artificial vertices";
        Arrangement testArrUpper;
        insert(testArrUpper, testUpper.begin(), testUpper.end());
        CGAL::draw(testArrLower);
        CGAL::draw(testArrUpper);
    }



    template<typename Kernel>
    typename K_visibility_region<Kernel>::Polygon K_visibility_region<Kernel>::find_visibility_region(int k, CGAL::Point_2<Kernel> p) {
        // transformations incase rotation is needed

        this->rotate = new Transformation(CGAL::ROTATION, sin(3.1415 / 100), cos(3.1415 / 100));
        if (p.y() == 0) {
            this->queryPoint = Point(p.x(), p.y() + 1);
            this->moveQueryPoint = new Transformation(CGAL::TRANSLATION, Vector(0, 1));
            std::cout << "moving query point" << std::endl;
            this->polygon = CGAL::transform(*this->moveQueryPoint, this->polygon);
        }
        else {
            this->queryPoint = Point(p.x(), p.y());
        }
        this->translateToOrigin = new Transformation(CGAL::TRANSLATION, Vector(-this->queryPoint.x(), -this->queryPoint.y()));
        this->translateBack = new Transformation(CGAL::TRANSLATION, Vector(this->queryPoint.x(), this->queryPoint.y()));


        auto x = this->polygon.vertex(0);
        std::cout << "query point: " << this->queryPoint << "\n";

        //  rotate polygon if horizontal line intersects a vertex
        while (this->isPointHorizontalWithVertex(this->queryPoint)) {
            std::cout << "Polygon has vertex on horizontal line through " << this->queryPoint << "\n";
            this->polygon = CGAL::transform(*translateToOrigin, this->polygon);
            this->polygon = CGAL::transform(*rotate, this->polygon);
            this->polygon = CGAL::transform(*translateBack, this->polygon);
            //auto test = (*this->translateBack) * (*this->rotate) * (*this->translateToOrigin);
            //this->polygon = CGAL::transform(test, this->polygon);
        }
        int i = 0;
        Polygon p2;
        for (auto e = this->polygon.edges_begin(); e != this->polygon.edges_end() - 1; ++e) {
            Segment t = *e;
            e[0].id() = e[0].id() == -1 ? i++ : e[0].id();
            e[1].id() = e[1].id() == -1 ? i++ : e[1].id();
            //    std::cout << "id: " << e->id() << "\n";
            i += 1;
        }
        /*i = 0;
        for (Point p : this->polygon.vertices()) {
            p.id() = i++;
            p2.push_back(p);
        }
        this->polygon = p2;*/
        std::vector<Segment> edges;
        for (Segment e : this->polygon.edges()) {
            edges.push_back(e);
        }

        for (Point point : this->polygon.vertices()) {
            //   Point tt(point.x(), point.y());
            std::cout << point.id() << std::endl;
        }




        assert(this->leftIntersectionIndex != -1);
        assert(this->rightIntersectionIndex != -1);

        getLowerUpper();
        getRadial(this->queryPoint);
       edges.insert(edges.end(), this->radialList.begin(), this->radialList.end());
        CGAL::draw(this->polygon);
        Arrangement test;
        insert(test, edges.begin(), edges.end());
      //  insert(test, Segment(Point(this->polygon.left_vertex()->x()- 1, this->queryPoint.y()), Point(this->polygon.right_vertex()->x() + 1, this->queryPoint.y())));
        // insert(arr, Segment(Point(0, 0)))
     //   Vertex_handle v = arr.insert_in_face_interior(Point(polygon.left_vertex()->x(), polygon.top_vertex()->y()), arr.unbounded_face());
        CGAL::draw(test);
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

        Point lm = *this->polygon.left_vertex();
        Point rm = *this->polygon.right_vertex();

        // Line line(p, Point(p.x() + 1, p.y()));
        Segment line(Point(lm.x() - 1, p.y()), Point(rm.x() + 1, p.y()));
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
    void K_visibility_region<Kernel>::radialHelper(Point p, Arrangement arr, int& artCounter) {
        Vert_decomp_list vd_list;
        CGAL::decompose(arr, std::back_inserter(vd_list));
        for (auto vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
            if (vd_iter->first->point().isArtificial()) {
                artCounter++;
                std::cout << "ARTIFICIAL" << std::endl;
            }
            Vertex_const_handle vertex = vd_iter->first;
            if (vertex->point().isArtificial()) {
                continue;
            }
            /*if (std::find(polygon.vertices_begin(), polygon.vertices_end(), inverseProjection(vertex->point())) == polygon.vertices_end()) {
                continue;
            }*/
            const Object_pair& curr = vd_iter->second;
            std::cout << "Vertex (" << vd_iter->first->point() << " w= " << vd_iter->first->point().hw() << ") -> (" << inverseProjection(vd_iter->first->point()) << " w= " << vd_iter->first->point().hw() << ")" << std::endl;
            Vertex_const_handle vh;
            Halfedge_const_handle hh;
            Face_const_handle fh;
            bool has_vh = false;
            bool has_hh = false;
            std::cout << " feature below: " << std::endl;

            if (!(has_hh = CGAL::assign(hh, curr.first)) && !(has_vh = CGAL::assign(vh, curr.first))) { // nothing below

                if (!(has_hh = CGAL::assign(hh, curr.second)) && !(has_vh = CGAL::assign(vh, curr.second))) { // nothing above
                    std::cout << "nothing above or below" << std::endl;
                    continue;
                }

                if (has_hh) { // nothing below, edge above
                    inv_projection_insert(hh, vertex);
                }
                else { // nothing below, vertex above
                    inv_projection_insert(vh, vertex);
                }
            }
            else { // something below
                bool has_hh2 = false;
                bool has_vh2 = false;

                Vertex_const_handle vh2;
                Halfedge_const_handle hh2;
                Face_const_handle fh2;
                if (!(has_hh2 = CGAL::assign(hh2, curr.second)) && !(has_vh2 = CGAL::assign(vh2, curr.second))) { // nothing above
                    std::cout << "nothing above" << std::endl;
                    if (has_hh) { // nothing above, edge below
                        inv_projection_insert(hh, vertex);
                    }
                    else { // nothing above, vertex below
                        inv_projection_insert(vh, vertex);
                    }
                    continue;
                }
                if (has_hh && has_hh2) {
                    inv_projection_insert(hh, hh2, vertex);
                }
                else if (has_vh && has_vh2) {
                    inv_projection_insert(vh, vh2, vertex);
                }
                else if (has_hh && has_vh2) {
                    inv_projection_insert(hh, vh2, vertex);
                }
                else if (has_vh && has_hh2) {
                    inv_projection_insert(vh, hh2, vertex);
                }
            }
        }
    }

    template <class Kernel>
    void K_visibility_region<Kernel>::getRadial(Point p) {
        
        int artCounter = 0;
        CGAL::insert(lowerArr, lowerEdgeList.begin(), lowerEdgeList.end());
        std::cout << "size of lowerArr: " << lowerEdgeList.size() << "\n";

        CGAL::draw(lowerArr);

        CGAL::insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
        CGAL::draw(upperArr);
        radialHelper(p, lowerArr, artCounter);
        radialHelper(p, upperArr, artCounter);
    }

    template <class Kernel>
    void K_visibility_region<Kernel>::projection_insert(std::vector<Segment>& vec, const Segment& seg) {
        Point s = this->projection(seg.start());
        Point e = this->projection(seg.end());
        Segment p(s, e);
        std::cout << "Segment: " << seg << " w = " << s.hw() << " " << e.hw() << ", projected: " << p << " w = " << s.hw() << " " << e.hw() << std::endl;
        vec.push_back(Segment(s, e));
    }

    template <class Kernel>
    void K_visibility_region<Kernel>::inv_projection_insert(const Segment& seg) {
        Point s = this->inverseProjection(seg.start());
        Point e = this->inverseProjection(seg.end());
        Segment p(s, e);
        std::cout << "Segment: " << seg << " w = " << s.hw() << " " << e.hw() << ", inv_projected: " << p << " w = " << s.hw() << " " << e.hw() << std::endl;
        this->radialList.push_back(Segment(s, e));
    }

   


    template <class Kernel>
    void  K_visibility_region<Kernel>::inv_projection_insert(Halfedge_const_handle &bhh, Halfedge_const_handle &uhh, Vertex_const_handle &qvh) {
        Point p1 = interpolate(bhh, qvh);
        Point p2 = interpolate(uhh, qvh);
        Segment s(p1, p2);
        this->radialList.push_back(s);
    }

    template <class Kernel>
    void  K_visibility_region<Kernel>::inv_projection_insert(Vertex_const_handle& bvh, Vertex_const_handle &uvh, Vertex_const_handle &qvh) {
        Segment s(bvh->point(), uvh->point());
        inv_projection_insert(s);
    }

    template <class Kernel>
    void  K_visibility_region<Kernel>::inv_projection_insert(Halfedge_const_handle &bhh, Vertex_const_handle& uvh, Vertex_const_handle &qvh) {
        Point p1 = interpolate(bhh, qvh);
        Point p2 = inverseProjection(qvh->point());
        Segment s(p1, p2);
        radialList.push_back(s);
    }

    template <class Kernel>
    void  K_visibility_region<Kernel>::inv_projection_insert(Vertex_const_handle &bvh, Halfedge_const_handle &uhh, Vertex_const_handle &qvh) {
        inv_projection_insert(uhh, bvh);
    }

    template <class Kernel>
    void  K_visibility_region<Kernel>::inv_projection_insert(Vertex_const_handle &bhh, Vertex_const_handle& qvh) {
        Segment s(bhh->point(), qvh->point());
        inv_projection_insert(s);
    }

    template <class Kernel>
    void  K_visibility_region<Kernel>::inv_projection_insert(Halfedge_const_handle &uhh, Vertex_const_handle &qvh) {
        Point p1 = interpolate(uhh, qvh);
        Point p2 = inverseProjection(qvh->point());
        Segment s(p1, p2);
        radialList.push_back(s);
    }

