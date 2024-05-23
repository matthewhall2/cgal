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
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Arr_naive_point_location.h>



#define EPSILON (1.0/1000.0)



    template<class Kernel>
    class K_visibility_region {

        using FT = typename Kernel::FT;
        using RT = typename Kernel::RT;
        using MK = typename K_Vis_Kernel<FT>;

        using Traits = typename CGAL::Arr_segment_traits_2<MK>;
        using GPS_Traits = typename CGAL::Gps_segment_traits_2<MK>;
        using Arrangement = CGAL::Arrangement_2<Traits>;
        using Point3 = CGAL::Point_3<Kernel>;
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
        using Face = typename Arrangement::Face;
        using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
        using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
        using Face_const_handle = typename Arrangement::Face_const_handle;
        using Face_iterator = typename Arrangement::Face_iterator;
        using X_monotone_curve = typename Arrangement::X_monotone_curve_2;
        using Face_handle = typename Arrangement::Face_handle;
        using Halfedge_iterator = typename Arrangement::Halfedge_iterator;
        using Edge_const_iterator = typename Arrangement::Edge_const_iterator;
        using Vertex_const_iterator = typename Arrangement::Vertex_const_iterator;

        using Halfedge_handle = typename Arrangement::Halfedge_handle;
        using RSPV = typename CGAL::Simple_polygon_visibility_2<Arrangement, CGAL::Tag_true>;



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

        Transformation* rotateForProj;

        Point queryPoint;
        Point lowestAboveL;
        Point highestBelowL;

        FT ymin, ymax, xmin, xmax;

        std::vector<std::pair<Point*, int>> intersectionPoints;
        int leftIntersectionIndex;
        int rightIntersectionIndex;
        Face_handle unboundedFace;

        void getLowerUpper();
        void lowerUpperHelper(int start, int end, int&, FT ly, FT uy, std::vector<Segment>&, std::vector<Segment>&);
        void merge();
        bool isPointHorizontalWithVertex(Point p);
        void findZeroVisibilty();
        void findNextRegion();
        int countFaceSides(Face_const_handle);
        void addHalfedge(Halfedge_const_handle stop);
        void insertBbox();
        void insertBboxHelper(Segment check, Segment stopCheck, int index);
        typedef enum {
           HORIZONTAL,
           VERTICAL
        } DIR;

        void findBbox(Point, DIR);
        Halfedge_handle extendRadial(Halfedge_const_handle);
        Vertex_handle splitEdge(Halfedge_handle, Point);

        Point projection(Point p) {
            // projection is the matrix
            /*1 0 0
            * 0 1 0
            * 0 1 -qp.y
            */
            //  std::cout << "projecting point: " << p << "\n";
       //     std::cout << "projecting with id: " << p.id() << std::endl;
            FT x = p.hx();
            FT y = p.hy();
            FT nz = CGAL::abs(y - this->queryPoint.y());
            //  std::cout << "query point: " << this->queryPoint << " hw = " << nz <<  std::endl;
            if (!p.isArtificial())             this->hws[p.id()] = nz;

            //    std::cout << "proj dif: " << nz << "\n";
            Point P(x / nz, y / nz, 1);
            //auto P = new Point(x , y , nz);
            P.isArtificial() = p.isArtificial();
            P.id() = p.id();
            assert(hws.find(-1) == hws.end());
            std::cout << "new w is " << nz << std::endl;
            std::cout << "y is " << y << std::endl;
            return  P;
        }

        Point inverseProjection(Point p) {
            std::cout << "inv projecting" << std::endl;
            // projection is the matrix
           /*1 0 0
           * 0 1 0
           * 0 1/qp.y -1/qp.y
           */
          //  std::cout << "inv proj id" << p.id() << std::endl;
         //   std::cout << "inv projecting point: " << p << "\n";
         //   std::cout << "homo: " << p.hx() << " " << p.hy() << std::endl;
            FT x = p.hx();
            FT y = p.hy();
            assert(p.id() != -1);
            assert(hws.find(p.id()) != hws.end());
            //  FT nz = y - this->queryPoint.y();
        //    FT nz = CGAL::inverse(this->queryPoint.y()) * y - CGAL::inverse(this->queryPoint.y()) * this->hws[p.id()];
          //  std::cout << "inv proj query point: " << this->queryPoint << " hw = " << this->hws[p.id()] << std::endl;

            //     std::cout << "proj dif: " << nz << "\n";
                 //std::cout << "query point: " << this->queryPoint << std::endl;
            std::cout << "homo " << hws[p.id()] << std::endl;
            Point P(x * this->hws[p.id()], y * this->hws[p.id()]);
            //auto P = new Point(p.x() * p.hw(), p.y() * p.hw(), 1);
            P.isArtificial() = p.isArtificial();
            assert(hws.find(-1) == hws.end());

            return  P;
        }

        Point3 p3(Point p) {
            return Point3(p.x(), p.y(), 0);
        }
        Polygon lowerProjected;
        Polygon upperProjected;
        void getRadial(Point p);
        void naiveRadial(Point p);
        void radialHelper(Point p, Arrangement arr, int &artCounter);
        
        //    CGAL::Unique_hash_map<std::tuple<FT, FT, FT, FT>, bool> ttt;
        //   boost::unordered_map<std::tuple<FT, FT, FT, FT>, std::vector<Point>> radialIntersectoinList;
          //  std::unordered_map<Point, bool> artificialVertexList;
        Arrangement upperArr;
        Arrangement lowerArr;
        Arrangement arr;
        Arrangement RegularVisibilityOutput;

        std::vector<Segment> upperEdgeList;
        std::vector<Segment> lowerEdgeList;
        std::vector<Segment> radialList;
        std::unordered_map<int, FT> hws;
        std::unordered_map<int, std::vector<Point>> intersectionList;
        Polygon testPolygon;

        void projection_insert(std::vector<Segment>& vec, const Segment& seg);

        
        Point interpolate(Halfedge_const_handle &hh, Vertex_const_handle &vh) {
            assert(!hh->is_fictitious());
            std::cout << "interpolating" << std::endl;
            Point source = hh->source()->point();
            Point target = hh->target()->point();
            Point p = vh->point();
            Line l = ((Segment)hh->curve()).supporting_line();
            Point pointOnLine(p.x(), l.y_at_x(p.x()));
            FT lineLength = CGAL::squared_distance(source, target);
            FT distFromSource = CGAL::squared_distance(source, pointOnLine);
            FT B = distFromSource / lineLength;
            B = CGAL::sqrt(B);
            assert(source.id() != -1);
            assert(target.id() != -1);
            std::cout << "ids " << source.id() << " " << target.id() << std::endl;
            assert(this->hws.find(source.id()) != hws.end());
            assert(this->hws.find(target.id()) != hws.end());
            auto s = hws.size();
            FT sourceW = this->hws[source.id()];
            assert(sourceW != 0);
            assert(s == hws.size());
            s = hws.size();
            FT targetW = this->hws[target.id()];
            FT a = (B / targetW) / (((1 - B) / sourceW) + (B / targetW));
            assert(targetW != 0);
            assert(s == hws.size());
            
            FT sx = source.x() * sourceW;
            FT sy = source.y() * sourceW;
            FT tx = target.x() * targetW;
            FT ty = target.y() * targetW;
            assert(polygon.has_on_boundary(Point(sx, sy)));
            assert(polygon.has_on_boundary(Point(tx, ty)));

            FT x = sx + a * (tx - sx);
            FT y = sy + a * (ty - sy);
            Point rp(x, y);
            assert(polygon.has_on_boundary(rp));
          //  Line l = (Line)hh->curve();
            return rp;
        }


        void inv_projection_insert(const Segment& seg);
        void inv_projection_insert(Halfedge_const_handle& bhh, Halfedge_const_handle& uhh, Vertex_const_handle &qvh);
        void inv_projection_insert(Vertex_const_handle &bvh, Vertex_const_handle &uvh, Vertex_const_handle &);

        void inv_projection_insert(Halfedge_const_handle &bhh, Vertex_const_handle &uhh, Vertex_const_handle &);
        void inv_projection_insert(Vertex_const_handle &bhh, Halfedge_const_handle &uhh, Vertex_const_handle &);

        void inv_projection_insert(Vertex_const_handle& bhh,Vertex_const_handle&);
        void inv_projection_insert(Halfedge_const_handle& uhh, Vertex_const_handle&);

        FT getAngle(Point p) {
            Point qp = this->queryPoint;
            FT angle = CGAL::approximate_angle(p3(Point(qp.x() + 1, qp.y())), p3(qp), p3(p));
            if ((p.x() < qp.x() && p.y() < qp.y()) || (p.x() > qp.x() && p.y() < qp.y())) {
                angle = -angle + 360;
            }
            return angle;
        }


        std::vector<std::pair<int, Halfedge_const_handle>> visibilityEdges; // id of full edge, and halfedge of visibilty region
        std::unordered_map<Halfedge_const_handle, int> testing;
        std::unordered_map<int, std::pair<Segment, Segment>> pointToSegments;

        typedef struct {
            bool hasMidpoint;
            Segment left;
            Segment right;
        } radialSegments;

        std::vector<Halfedge_const_handle> currentRegionHalfedges;
        std::vector<Halfedge_const_handle> nextRegionHalfedges;
        std::vector<Segment> bboxEdges;
        std::vector<Point> bboxCorners;

       
    };

    template<typename Kernel>
    K_visibility_region<Kernel>::K_visibility_region(CGAL::Polygon_2<Kernel> p) {
        this->polygon.clear();
        if (p.orientation() == CGAL::CLOCKWISE) {
            p.reverse_orientation();
        }
        this->polygon.clear();
        int i = 0;
        for (CGAL::Point_2<Kernel> point : p.vertices()) {
            Point np(point.x(), point.y());
            np.id() = i;
            CGAL::Segment_2<Kernel> e1 = p.edge(i);
            CGAL::Segment_2<Kernel> e2 = (i - 1 >= 0 ? p.edge(i - 1) : p.edge(p.edges().size() - 1));
            Kernel::FT res = CGAL::determinant(e2.to_vector(), e1.to_vector());
            if (res < 0) {
                std::cout << "reflex" << std::endl;
                np.isReflex() = true;
            }

            this->polygon.push_back(np);
            i += 1;
        }
        CGAL::draw(p);


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
     //   std::cout << "count in test: " << count << std::endl;


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

        Arrangement test;
        insert_non_intersecting_curve(test, Segment(Point(0, 0), Point(0, 1)));
        insert_non_intersecting_curve(test, Segment(Point(0, 1), Point(1, 1)));
        insert_non_intersecting_curve(test, Segment(Point(1, 1), Point(1, 0)));
        Halfedge_handle edge = insert_non_intersecting_curve(test, Segment(Point(1, 0), Point(0, 0)));
        Face_handle f = (edge->face()->is_unbounded() ? edge->twin()->face() : edge->face());
       // test.insert_in_face_interior(X_monotone_curve(Point(0.5, 0), Point(0.5, 1)), f);
     //   CGAL::draw(test);
      //  insert_non_intersecting_curve(test, Segment(Point(0.5, 0), Point(0.5, 1)));
       // insert(test, Segment(Point(0, 0), Point(1, 1)));
    }

    template<typename Kernel>
    void K_visibility_region<Kernel>::lowerUpperHelper(int start, int end, int &artCounter, FT ly, FT uy, std::vector<Segment> &testLower, std::vector<Segment> &testUpper) {
        static enum {
            ABOVE = 0,
            BELOW = 1
        } state = (start + 1 < polygon.vertices().size() ? (polygon.vertex(start).y() < polygon.vertex(start + 1).y() ? BELOW : ABOVE) : (polygon.vertex(start).y() < polygon.vertex(0).y() ? BELOW : ABOVE));
        Segment e;
        Line l;
        int nv = this->polygon.vertices().size();

        for (int i = start; i < end; i++) {
            e = this->polygon.edge(i);// *(this->polygon.edges_begin() + i);
            std::cout << "dividing: edge has ids " << e[0].id() << " " << e.target().id() << " artificial?: " << (e.source().isArtificial() ? "true" : "false") << " " << (e.target().isArtificial() ? "true" : "false") << std::endl;

            l = e.supporting_line();
            if (state == ABOVE) {
               // std::cout << this->polygon.edge(i) << "\n";
                if (this->intersectionPoints[i].first == nullptr) {
                    projection_insert(upperEdgeList, e);
                    testUpper.push_back(e);
                    continue;
                }
                Point low(l.x_at_y(ly), ly);
                low.id() = e.start().id();
                artCounter++;
                low.isArtificial() = true;
                Segment lowSeg(low, e.end());
                // its only inserting artificial vertices in the list? only articifial have valid ids?

                lowSeg.id() = e.id();
                projection_insert(lowerEdgeList, lowSeg);
                testLower.push_back(lowSeg);

                Point up(l.x_at_y(uy), uy);
                up.isArtificial() = true;
                up.id() = e.end().id();
                 artCounter++;
                Segment upSeg(e.start(), up);

                upSeg.id() = e.id();
                projection_insert(upperEdgeList, upSeg);
                testUpper.push_back(upSeg);

                state = BELOW;
            }
            else if (state == BELOW) {

                if (this->intersectionPoints[i].first != nullptr) {
                 //   std::cout << *this->intersectionPoints[i].first << "\n";

                    Point low(l.x_at_y(ly), ly);
                    low.isArtificial() = true;
                    low.id() = e.start().id();
                    artCounter++;
                    Segment lowSeg(e.start(), low);
                    lowSeg.id() = e.id();
                    projection_insert(lowerEdgeList, lowSeg);
                    testLower.push_back(lowSeg);

                    Point up(l.x_at_y(uy), uy);
                    up.isArtificial() = true;
                    up.id() = e.end().id();
                    artCounter++;

                    Segment upSeg(up, e.end());
                    upSeg.id() = e.id();

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
        int artCounter = 0;
        int nv = this->polygon.vertices().size();
     //   assert(this->queryPoint.y() - highestBelowL.y() > 0);
     //   assert(lowestAboveL.y() - this->queryPoint.y() > 0);

        FT ly = this->queryPoint.y() - (CGAL::abs(this->queryPoint.y() - highestBelowL.y())) * EPSILON;
        FT uy = this->queryPoint.y() + (CGAL::abs(this->lowestAboveL.y() - this->queryPoint.y())) * EPSILON;
   
        lowerEdgeList.clear();
        upperEdgeList.clear();
   
        this->lowerUpperHelper(this->leftIntersectionIndex, this->intersectionPoints.size(), artCounter, ly, uy, testLower, testUpper);
        this->lowerUpperHelper(0, this->leftIntersectionIndex, artCounter, ly, uy, testLower, testUpper);

     //   std::cout << "trgethsrh";
        std::cout << "size is " << lowerEdgeList.size() << std::endl;
        std::cout << "number of homogeneous: " << this->hws.size() << std::endl;
        Arrangement testArrLower;
        insert(testArrLower, testLower.begin(), testLower.end());
        std::cout << "There are " << artCounter << " artificial vertices";
        Arrangement testArrUpper;
        insert(testArrUpper, testUpper.begin(), testUpper.end());
      //  CGAL::draw(testArrLower);
     //   CGAL::draw(testArrUpper);
    }

    template<typename Kernel>
    typename K_visibility_region<Kernel>::Polygon K_visibility_region<Kernel>::find_visibility_region(int k, CGAL::Point_2<Kernel> p) {
        if (k < 0) {
            std::stringstream ss;
            ss << "k must be non-negative";
            throw std::invalid_argument(ss.str());
        }
       

        CGAL_precondition(k >= 0);
        visibilityEdges.clear();
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
        rotateForProj = new Transformation(CGAL::ROTATION, 1, 0);

        auto x = this->polygon.vertex(0);
        std::cout << "query point: " << this->queryPoint << "\n";

        //  rotate polygon if horizontal line intersects a vertex
        while (this->isPointHorizontalWithVertex(this->queryPoint)) {
            std::cout << "Polygon has vertex on horizontal line through " << this->queryPoint << "\n";
            this->polygon = CGAL::transform(*translateToOrigin, this->polygon);
            this->polygon = CGAL::transform(*rotate, this->polygon);
            this->polygon = CGAL::transform(*translateBack, this->polygon);
        }
        assert(this->leftIntersectionIndex != -1);
        assert(this->rightIntersectionIndex != -1);

        for (int i = 0; i < polygon.vertices().size(); i++) {
            pointToSegments[polygon.vertex(i).id()] = std::make_pair(polygon.edge(i), polygon.edge(i));
        }
        CGAL::draw(polygon);
        getLowerUpper();
       // getRadial(this->queryPoint);
        insert_non_intersecting_curves(arr, polygon.edges().begin(), polygon.edges().end());
        insertBbox();
        naiveRadial(this->queryPoint);
     //   insert_non_intersecting_curves(arr, radialList.begin(), radialList.end());
        findZeroVisibilty();

        Polygon region;
        for (Halfedge_const_handle h : nextRegionHalfedges) {
            region.push_back(h->source()->point());
        }
        CGAL::draw(region);
        while (k > 0) {
            std::cout << "finding " << k << "-visibility region" << std::endl;
            findNextRegion();
            region.clear();
            for (Halfedge_const_handle h : nextRegionHalfedges) {
                region.push_back(h->source()->point());
            }
            CGAL::draw(region);
            k -= 1;
        }

      

      /*  CGAL::draw(this->polygon);
        CGAL::draw(RegularVisibilityOutput);
        Arrangement testarr;*/
       
       // CGAL::draw(arr);
     //   Vertex_handle v = arr.insert_in_face_interior(Point(polygon.left_vertex()->x(), polygon.top_vertex()->y()), arr.unbounded_face());
        int test = 5;
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
          
            const Object_pair& curr = vd_iter->second;
            Vertex_const_handle vh;
            Halfedge_const_handle hh;
            Face_const_handle fh;
            bool has_vh = false;
            bool has_hh = false;

            std::cout << " feature below: " << std::endl;
            //none intersect a vertex in current example
            if (!(has_hh = CGAL::assign(hh, curr.first)) && !(has_vh = CGAL::assign(vh, curr.first))) { // nothing below

                if (!(has_hh = CGAL::assign(hh, curr.second)) && !(has_vh = CGAL::assign(vh, curr.second))) { // nothing above
                    std::cout << "nothing above or below" << std::endl;
                    continue;
                }

                if (has_hh) { // nothing below, edge above
                    inv_projection_insert(hh, vertex);
                }
                else { // nothing below, vertex above
                  //  inv_projection_insert(vh, vertex);
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
                    //    inv_projection_insert(vh, vertex);
                    }
                    continue;
                }
                if (has_hh && has_hh2) {
                    inv_projection_insert(hh, hh2, vertex);
                }
                else if (has_vh && has_vh2) {
                  //  inv_projection_insert(vh, vh2, vertex);
                }
                else if (has_hh && has_vh2) {
                   // inv_projection_insert(hh, vh2, vertex);
                }
                else if (has_vh && has_hh2) {
                   // inv_projection_insert(vh, hh2, vertex);
                }
            }
        }
    }

    template <class Kernel>
    void K_visibility_region<Kernel>::getRadial(Point p) {
        testPolygon.clear();
        int artCounter = 0;
        CGAL::insert(lowerArr, lowerEdgeList.begin(), lowerEdgeList.end());
        std::cout << "size of lowerArr: " << lowerEdgeList.size() << "\n";

        CGAL::draw(lowerArr);

        CGAL::insert(upperArr, upperEdgeList.begin(), upperEdgeList.end());
       CGAL::draw(upperArr);

        std::vector<Segment> el2;
        std::vector<Segment> ul2;
        for (Segment s : lowerEdgeList) {
            if (s.start().isArtificial() || s.end().isArtificial()) {
                continue;
            }
            el2.push_back(s);
        }
        for (Segment s : upperEdgeList) {
            if (s.start().isArtificial() || s.end().isArtificial()) {
                continue;
            }
            ul2.push_back(s);
        }

        Arrangement arr1;
        Arrangement arr2;
        CGAL::insert(arr1, el2.begin(), el2.end());
    //    CGAL::draw(arr1);
        CGAL::insert(arr2, ul2.begin(), ul2.end());
    //    CGAL::draw(arr2);


        radialHelper(p, lowerArr, artCounter);
        radialHelper(p, upperArr, artCounter);
        std::cout << "done" << std::endl;
      //  CGAL::draw(testPolygon);
    }

    template <class Kernel>
    void K_visibility_region<Kernel>::projection_insert(std::vector<Segment>& vec, const Segment& seg) {
        Point s = this->projection(seg.start());
        Point e = this->projection(seg.end());
        Segment p(s, e);
    
        vec.push_back(p);
    }

    template <class Kernel>
    void K_visibility_region<Kernel>::inv_projection_insert(const Segment& seg) {
        Point s = this->inverseProjection(seg.start());
        Point e = this->inverseProjection(seg.end());
        Segment p(s, e);
        this->radialList.push_back(p);
    }

    template <class Kernel>
    void  K_visibility_region<Kernel>::inv_projection_insert(Halfedge_const_handle &bhh, Halfedge_const_handle &uhh, Vertex_const_handle &qvh) {
        Point p1 = interpolate(bhh, qvh);
        Point p2 = interpolate(uhh, qvh);
        Segment s(p1, p2);
       // std::cout << s << std::endl;
       assert(!s.is_degenerate());
      
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

    template <class Kernel>
    void K_visibility_region<Kernel>::naiveRadial(Point p) {
#define sid(seg) (seg.start().id() + seg.end().id())
        Point left, right;
        Segment leftSeg, rightSeg;
        FT leftDist, rightDist;
        int pid = polygon.vertices().size();
        for (auto v = arr.vertices_begin(); v != arr.vertices_end(); v++) {
            Point vertex = v->point();
            leftDist = FT(1000000000);
            rightDist = FT(100000000);
            bool foundLeft = false;
            bool foundRight = false;
            bool isVertical = vertex.x() == p.x();
            bool side = isVertical ? vertex.y() < p.y() : vertex.x() < p.x();
            Line ray(vertex, p);
            for (auto e = arr.edges_begin(); e != arr.edges_end(); e++) {
                Segment edge = e->curve();
                Point intPoint;
                const auto result = CGAL::intersection(ray, edge);
                if (!result) continue;

                if (const Point* ip = boost::get<Point>(&*result)) {
                    if (*ip == v->point()) {
                        continue;
                    }

                    intPoint = *ip;
                }
                else {
                    continue;
                    /*const Segment* s = boost::get<Segment>(&*result);
                    FT d1 = CGAL::squared_distance(vertex, s->start());
                    FT d2 = CGAL::squared_distance(vertex, s->end());
                    
                    intPoint = d1 < d2 ? s->start() : s->end();*/
                }

                bool side2 = isVertical ? intPoint.y() < p.y() : intPoint.x() < p.x();
                if (side != side2) continue;

                if (isVertical ? intPoint.y() < vertex.y() : intPoint.x() < vertex.x()) { // left side
                    foundLeft = true;
                    FT dist = CGAL::squared_distance(vertex, intPoint);
                    if (dist < leftDist) {
                        left = intPoint;
                        leftDist = dist;
                        leftSeg = edge;
                    }
                    
                }
                else if (isVertical ? intPoint.y() > vertex.y() : intPoint.x() > vertex.x()) {
                    foundRight = true;
                    FT dist = CGAL::squared_distance(vertex, intPoint);
                    if (dist < rightDist) {
                        right = intPoint;
                        rightDist = dist;
                        rightSeg = edge;
                    }
                }
            }
            // create structure to store segments and midpoint if applicatable
            // when travering, check if there is a midterm (which means had relfex vertex)
            // can use intert_in_face_interior to speed up (already have face from edge->face())
            // also have map from new point id (new point being point in radial) to the segment id if itersects
            // if vertex, have 2 ids
            // if both point ids are in range (0, size of polygon), then both are vertices
            // if either point is a vertex, check if other point is on same segment
            // if yes, has intersections
            // if not, no possible intersections

            // if vertex, face is already there, so loop and add halfedges until target is target

            //use face.isunbounded to check is on outer halfedge
            

            if (foundLeft && foundRight) {
                if (vertex.isReflex()) {
                    left.isTip() = true;
                    right.isTip() = true;
                }
                left.baseAwayFromQP() = side;
                right.baseAwayFromQP() = !side;
                right.id() = pid++;
                left.id() = pid++;

                intersectionList[sid(rightSeg)].push_back(right);
                intersectionList[sid(leftSeg)].push_back(left);
               // insert(arr, Segment(left, vertex));
                radialList.push_back(Segment(left, right));
              //  CGAL::draw(arr);
              //  insert_non_intersecting_curve(arr, Segment(right, vertex));

                // split vertex to be able to insert using insert_non_intersecting_curves()
               // radialList.push_back(Segment(left, vertex));
              //  radialList.push_back(Segment(vertex, right));
                continue;
            }
            


            if (foundRight) {
                if (vertex.isReflex()) {
                    right.isTip() = true;
                }
                right.baseAwayFromQP() = !side;
                right.id() = pid++;

                radialList.push_back(Segment(vertex, right));
                intersectionList[sid(rightSeg)].push_back(right);
              /*  if (isVertical ? vertex.y() > p.y() : vertex.x() > p.x()) {
                    radialList.push_back(Segment(p, vertex));
                }*/
                continue;
                
            }

            if (foundLeft) {
                if (vertex.isReflex()) {
                    left.isTip() = true;
                }
                left.baseAwayFromQP() = side;
                left.id() = pid++;

                radialList.push_back(Segment(left, vertex));
                intersectionList[sid(leftSeg)].push_back(left);
               /* if (isVertical ? vertex.y() < p.y() : vertex.x() < p.x()) {
                    radialList.push_back(Segment(p, vertex));
                }*/
                continue;
            }

            /*if (isVertical ? vertex.y() > p.y() : vertex.x() > p.x()) {
                radialList.push_back(Segment(p, vertex));
            }
            if (isVertical ? vertex.y() < p.y() : vertex.x() < p.x()) {
                radialList.push_back(Segment(p, vertex));
            }*/
            


        }
    }

    template<class Kernel>
    void K_visibility_region<Kernel>::findZeroVisibilty() {
      //  CGAL::draw(arr);
      
        CGAL::draw(arr);
        Halfedge_iterator h = arr.halfedges_begin();
        Face_handle f = h->face()->is_unbounded() ? h->twin()->face() : h->face();


        Arrangement polyArr;
        insert(polyArr, polygon.edges_begin(), polygon.edges_end());
        Face_const_handle inputPolygonFace, radialDecompZeroFace;
        CGAL::Arr_naive_point_location<Arrangement> inputPolyPointLocation(polyArr);
        CGAL::Arr_point_location_result<Arrangement>::Type inputFaceResult = inputPolyPointLocation.locate(queryPoint);
        // The query point locates in the interior of a face
        inputPolygonFace = *boost::get<Face_const_handle>(&inputFaceResult);
        RegularVisibilityOutput.clear();
        RSPV regular_visibility(polyArr);
        regular_visibility.compute_visibility(queryPoint, inputPolygonFace, RegularVisibilityOutput);
        
        CGAL::Arr_naive_point_location<Arrangement> zeroFacePointLocation(RegularVisibilityOutput);
        CGAL::Arr_point_location_result<Arrangement>::Type zeroFaceResult = zeroFacePointLocation.locate(queryPoint);
        Face_const_handle zeroVisibiltyFace = *boost::get<Face_const_handle>(&zeroFaceResult);

        insert(arr, radialList.begin(), radialList.end());
        CGAL::Arr_naive_point_location<Arrangement> radialArrPointLocation(arr);
        CGAL::Arr_point_location_result<Arrangement>::Type radialDecompZeroResult = radialArrPointLocation.locate(queryPoint);
        radialDecompZeroFace = *boost::get<Face_const_handle>(&radialDecompZeroResult);
        assert(*radialDecompZeroFace == *f);

        auto edge = radialDecompZeroFace->outer_ccb();
        Halfedge_const_handle edge2 = edge->twin();
        assert(radialDecompZeroFace == edge->face());
        std::cout << typeid(Point(1, 1)).name() << std::endl;
        std::cout << typeid(*edge).name() << std::endl;
        Face_handle fh = &(*(edge->face()));
        for (auto testFace = arr.faces_begin(); testFace != arr.faces_end(); testFace++) {
            if (&*testFace == &*radialDecompZeroFace) {
                std::cout << "found match if face" << std::endl;
            }
        }
        polyArr.insert_in_face_interior(Point(1, 1), fh);
        assert(&*fh == &*(edge->face()));
       // edge = edge2;
       // assert(face2->inner_ccbs_begin() != face2->inner_ccbs_end());
        do {
            std::cout << "reg vis : [" << edge->source()->point() << " -> " << edge->target()->point() << "]" << std::endl;
            int sid = edge->source()->point().id();
            int tid = edge->target()->point().id();
            Segment curr = edge->curve();
            Segment nex = edge->next()->curve();
            if (curr.direction() == nex.direction()) {
                std::cout << "same dir as next edge" << std::endl;
            }
            assert(sid != -1 && tid != -1);
            assert(radialDecompZeroFace == edge->face());
         //   FT a = getAngle(edge->source()->point());
            Segment e = (Segment)edge->curve();
            Point p = e.start();
            testing[edge] = 1;
            nextRegionHalfedges.push_back(edge);
          //  visibilityEdges.push_back(edge->source()->point());
         //   std::cout << "angle is: " << a << std::endl;
            std::cout << "id is " << edge->source()->point().id() << std::endl;
            if (edge->source()->point().id() == -1) {
                Halfedge_const_handle hh;
                CGAL::Arr_point_location_result<Arrangement>::Type obj2 = inputPolyPointLocation.locate(edge->source()->point());
                hh = *boost::get<Halfedge_const_handle>(&obj2);
            }
        //  std::cout << ((Segment)(edge->curve())).start().id() << std::endl;
          //edge = edge->next();
            edge++;
        } while (edge != radialDecompZeroFace->outer_ccb());
       // std::reverse(visibilityEdges.begin(), visibilityEdges.end());
    }

    template<class Kernel>
    void K_visibility_region<Kernel>::findNextRegion() {
        auto s = nextRegionHalfedges.size();
        currentRegionHalfedges = nextRegionHalfedges;
        nextRegionHalfedges.clear();
        assert(currentRegionHalfedges.size() == s);
        assert(nextRegionHalfedges.size() == 0);
        Arrangement test;
        for (Halfedge_const_handle edge : currentRegionHalfedges) {
            test.clear();
            Halfedge_const_handle oppFaceEdge = edge->twin();
            //   Halfedge_const_handle queryEdge = oppFaceEdge->next();
            Face_handle face = &*(oppFaceEdge->face());
            //    Segment oppFaceEdgeCurve = oppFaceEdge->curve();
            //    Segment nextCurve = queryEdge->curve();
                // will be different for subregion of p (not plane)
            if (face->is_unbounded()) {
                nextRegionHalfedges.push_back(edge);
                // addHalfedge(edge, oppFaceEdge);
                // nextRegionHalfedges.push_back(edge);
                continue;
            }
            // if not tip, no new edge to insert, just have to add halfedge to list
            // if tip, the base could be 
           // if (!oppFaceEdge->target()->point().isTip() && !oppFaceEdge->source()->point().isTip()) { // halfedge already in arrangement
            addHalfedge(oppFaceEdge);
            //nextRegionHalfedges.push_back(queryEdge);
      //  }
      //  else { // traverse face (at most 4 edges) to find intersection

         //   extendRadial(edge);

        //}
            test.clear();
            for (Halfedge_const_handle h : nextRegionHalfedges) {
                insert(test, Segment(h->source()->point(), h->target()->point()));
            }
        //    CGAL::draw(test);
        }
    }

    template<class Kernel>
    void K_visibility_region<Kernel>::addHalfedge(Halfedge_const_handle edge) {
        Vertex_handle target = &*edge->target();
        Halfedge_const_handle toInsert;
        bool isTip = target->point().isTip();
        if (isTip == false) {
            toInsert = edge->next();
        }
        else {
            bool lidAway = target->point().baseAwayFromQP();
            if (lidAway == false) {
                toInsert = extendRadial(edge);
            }
            else {
                Segment s1 = edge->next()->curve();
                Segment s2 = edge->next()->next()->curve();

                toInsert = (s1.direction() == s2.direction() ? edge->next()->next() : edge->next());
            }
        }

        if (nextRegionHalfedges.size() == 0) {
            nextRegionHalfedges.push_back(toInsert);
        }
        else {
            Halfedge_const_handle last = nextRegionHalfedges.back();
            //   std::cout << "last edge target " << last->target()->point() << ", edge target " << edge->target()->point() << std::endl;
            if (toInsert->target()->point() != last->target()->point()) {
                nextRegionHalfedges.push_back(toInsert);
            }
            else {
                toInsert = toInsert->next();
                nextRegionHalfedges.push_back(toInsert);
            }
            // not called when tip, so not needed
            /*else if (((Segment)(edge->curve())).direction() == ((Segment)(edge->next()->curve())).direction()) {
                nextRegionHalfedges.push_back(edge->next());
                edge = edge->next();
            }*/
        }
        int count = 0;
        //(next->source()->point().isReflex() && !next->target()->point().baseAwayFromQP())
        Halfedge_const_handle next = toInsert->next();
        while (((next->target()->point() != edge->source()->point()) && (!next->source()->point().isReflex()))) {
            nextRegionHalfedges.push_back(next);
            next = next->next();
            count += 1;
        }
        std::cout << "added " << count << " edges" << std::endl;

    }

    template<class Kernel>
    int K_visibility_region<Kernel>::countFaceSides(Face_const_handle f) {
        int count = 0;
        auto cc = f->outer_ccb();
        do {
            count++;
            cc++;
        } while (cc != f->outer_ccb());
        return count;
    }


    //only target matters?
    // if base away from qp, add edge->next()->next()
    // unless edge->next->face != 
    // stop when target of newly inserted (into regionnlisy) edge is either reflex, or source of stop point
    template<class Kernel>
    typename K_visibility_region<Kernel>::Halfedge_handle  K_visibility_region<Kernel>::extendRadial(Halfedge_const_handle edge) {
        Point target = edge->target()->point();
        std::cout << "tip" << std::endl;
        Halfedge_const_handle queryEdge = edge->next();
        Line extender(queryPoint, target);
        Halfedge_handle r;
        Segment nextCurve;
        while (queryEdge != edge) {
            std::cout << "test";
            //  assert(queryEdge->face() == face);
            nextCurve = queryEdge->curve();
            const auto result = CGAL::intersection(extender, nextCurve);
            if (!result) {
                queryEdge = queryEdge->next();

                continue;
            }

            if (const Point* p = boost::get<Point>(&*result)) {
                
                if (*p == queryEdge->source()->point()) {
                    // queryEdge = queryEdge->next();
                    return &*queryEdge->prev();
                }
                if (*p == queryEdge->target()->point()) {
                    return &*queryEdge->next();
                }
                Vertex_handle split = splitEdge(&*queryEdge, *p);
                Halfedge_handle newEdge = arr.insert_at_vertices(X_monotone_curve(edge->target()->point(), split->point()), (Vertex_handle)(&*edge->target()), split);
                //   Halfedge_const_handle newEdge = arr.insert_in_face_interior(X_monotone_curve(start, *p), face);
                Segment newCurve = newEdge->curve();
                if (newCurve.direction() == extender.direction()) {
                    r = newEdge;
                    //  nextRegionHalfedges.push_back(newEdge);
                }
                else {
                    // nextRegionHalfedges.push_back(newEdge->twin());
                    r = newEdge->twin();
                }
                break;
            }
            queryEdge = queryEdge->next();
        }
        return r;
    }

    template<class Kernel>
    typename K_visibility_region<Kernel>::Vertex_handle K_visibility_region<Kernel>::splitEdge(Halfedge_handle edge, Point p) {
        p.isTip() = true;
        Segment c1(edge->source()->point(), p);
        Segment c2(p, edge->target()->point());
        Halfedge_handle split = arr.split_edge(edge, c1, c2);
        return split->target();
    }

    template<class Kernel>
    void K_visibility_region<Kernel>::insertBbox() {
        if (polygon.orientation() == CGAL::CLOCKWISE) {
            polygon.reverse_orientation();
        }
        bboxEdges.clear();
        Point left = *polygon.left_vertex();
        Point right = *polygon.right_vertex();
        Point top = *polygon.top_vertex();
        Point bottom = *polygon.bottom_vertex();

        this->ymin = bottom.y(), this->xmin = left.x(), this->ymax = top.y(), this->xmax = right.y();
        int i = 0;
        int leftIndex, rightIndex, topIndex, bottomIndex;
        for (Point p : polygon.vertices()) {
            if (p == left) {
                leftIndex = i;
            }
            else if (p == right) {
                rightIndex = i;
            }
            else if (p == top) {
                topIndex = i;
            }
            else if (p == bottom) {
                bottomIndex = i;
            }
            i++;
        }

        Point topRight(right.x(), top.y());
        Point topLeft(left.x(), top.y());
        Point bottomLeft(left.x(), bottom.y());
        Point bottomRight(right.x(), bottom.y());
        bboxCorners.push_back(topRight);
        bboxCorners.push_back(topLeft);
        bboxCorners.push_back(bottomRight);
        bboxCorners.push_back(bottomRight);

        Segment topSeg(topRight, topLeft);
        Segment leftSeg(topLeft, bottomLeft);

        Segment bottomSeg(bottomLeft, bottomRight);
        Segment rightSeg(bottomRight, topRight);

        std::vector<Segment> edges;

        insertBboxHelper(topSeg, leftSeg, topIndex);
        insertBboxHelper(leftSeg, bottomSeg, leftIndex);
        insertBboxHelper(bottomSeg, rightSeg, bottomIndex);
        insertBboxHelper(rightSeg, topSeg, rightIndex);

        for (Segment s : bboxEdges) {
            assert(s.is_horizontal() || s.is_vertical());
        }
        insert(arr, bboxEdges.begin(), bboxEdges.end());
    }

    template<class Kernel>
    void K_visibility_region<Kernel>::insertBboxHelper(Segment checkedEdge, Segment stopCheck, int index) {
        Point lastPoint = checkedEdge.source();
        auto size = polygon.vertices().size();
        Segment edge = polygon.edge(index);
        auto  result = CGAL::intersection(edge, stopCheck);
        while (!result) {
            auto const intPoint = CGAL::intersection(edge, checkedEdge);
            if (!intPoint) {
                index += 1;
                if (index >= size) {
                    index = 0;
                }
                edge = polygon.edge(index);
                result = CGAL::intersection(edge, stopCheck);
                continue;
            }

            if (const Point* p = boost::get<Point>(&*intPoint)) {
                if (*p != lastPoint) {
                    bboxEdges.push_back(Segment(lastPoint, *p));
                    lastPoint = *p;
                }
            }
            else if (const Segment* s = boost::get<Segment>(&*intPoint)) {
                if (s->source() != lastPoint) {
                    bboxEdges.push_back(Segment(lastPoint, s->source()));
                    lastPoint = s->target();

                }
            }
            index += 1;
            if (index >= size) {
                index = 0;
            }
            edge = polygon.edge(index);
            result = CGAL::intersection(edge, stopCheck);

        }
        if (bboxEdges.back().target() != checkedEdge.target()) {
            bboxEdges.push_back(Segment(lastPoint, checkedEdge.target()));
        }
    }

    template<class Kernel>
    void K_visibility_region<Kernel>::findBbox(Point p, DIR d) {
        Segment s1 = DIR == HORIZONTAL ? Segment(this->xmin, p.y()) : Segment(p.x(), this->ymin);
        Segment s2 = DIR == HORIZONTAL ? Segment(this->xmax, p.y()) : Segment(p.x(), this->ymax);


    }



