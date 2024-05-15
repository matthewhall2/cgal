// Author(s) : Camille Wormser, Pierre Alliez

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>

typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;

typedef std::list<Triangle>::const_iterator Iterator;
typedef CGAL::AABB_triangle_primitive_3<K, Iterator> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

int main()
{
  Point a(1.0, 0.0, 0.0);
  Point b(0.0, 1.0, 0.0);
  Point c(0.0, 0.0, 1.0);
  Point d(0.0, 0.0, 0.0);

  std::list<Triangle> triangles;
  triangles.push_back(Triangle(a, b, c));
  triangles.push_back(Triangle(a, b, d));
  triangles.push_back(Triangle(a, d, c));

  // constructs AABB tree
  Tree tree(triangles.begin(), triangles.end());

  // counts #intersections
  Ray ray_query(a, b);
  assert(tree.number_of_intersected_primitives(ray_query) == 3);

  // compute closest point and squared distance
  Point point_query(3.0, 2.0, 2.0);
  Point closest_point = tree.closest_point(point_query);
  assert(closest_point == a);
  assert(tree.squared_distance(point_query) == 12);

  return EXIT_SUCCESS;
}
