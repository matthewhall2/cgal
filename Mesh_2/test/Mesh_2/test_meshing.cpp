// 154 515 565
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/IO/File_poly.h>

#include <fstream>
#include <iostream>

template <class CTr>
typename CTr::size_type number_of_constrained_edges(const CTr& tr)
{
  typename CTr::size_type nedges = 0;
  for(typename CTr::Finite_edges_iterator eit = tr.finite_edges_begin();
      eit != tr.finite_edges_end();
      ++eit)
    if(tr.is_constrained(*eit))
      ++nedges;
  return nedges;
}

template <typename K>
struct Tester {
  void operator()() const {

    typedef CGAL::Triangulation_vertex_base_2<K> Vb;
    typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
    typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
    typedef typename CDT::size_type size_type;

    typedef typename CDT::Point Point;

    CDT cdt;

    std::vector<Point> seeds;
    seeds.reserve(32);

    std::cerr << "Reading fish-and-rectangle.poly...";
    std::ifstream poly_file("fish-and-rectangle.poly");
    CGAL::read_triangle_poly_file(cdt, poly_file, std::back_inserter(seeds));

    const size_type inititial_number_of_vertices = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of seeds: " << seeds.size() << "\n\n";

    std::cerr << "Saving the triangulation...\n\n";
    CDT cdt2 = cdt;

    std::cerr << "1/ First tests:\n\n";

    std::cerr << "Meshing the triangulation with size 0...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria());
    const size_type number_of_vertices0 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    CGAL_assertion( 64 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 72 );
    CGAL_assertion( seeds.size() == 3 );

    std::cerr << "Meshing the triangulation with size 0.2...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria(0.125, 0.2));

    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    CGAL_assertion( 190 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 210 );

    std::cerr << "Meshing the triangulation with size 0.1...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria(0.125, 0.1));
    const size_type number_of_vertices1 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
    CGAL_assertion( 580 <= cdt.number_of_vertices() &&
                    cdt.number_of_vertices() <= 640 );

    cdt = cdt2;
    std::cerr << "Triangulation restored.\n";
    std::cerr << "Number of vertices: " << cdt.number_of_vertices() << "\n\n";

    std::cerr << "Meshing the triangulation with Delaunay_mesh_criteria_2<CDT>()...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 CGAL::Delaunay_mesh_criteria_2<CDT>());
    const size_type number_of_vertices0bis = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() << "\n\n";
  
    CGAL_assertion( number_of_vertices0 == number_of_vertices0bis );

    cdt = cdt2;
    std::cerr << "Triangulation restored.\n";
    std::cerr << "Number of vertices: " << cdt.number_of_vertices() << "\n\n";

    std::cerr << "2/ Comparaison between refine_Delaunay_mesh_2() and other"
              << " possibilities:\n\n";

    std::cerr << "Meshing the triangulation with size 0.1, with "
              << "refine_Delaunay_mesh_2()...";
    CGAL::refine_Delaunay_mesh_2(cdt,
                                 seeds.begin(), seeds.end(),
                                 Criteria(0.125, 0.1));
    const size_type number_of_vertices1bis = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\n\n";

    CGAL_assertion( number_of_vertices1bis <= number_of_vertices1 );

    cdt = cdt2;

    std::cerr << "Meshing the triangulation with size 0.1, with "
              << "mesher.refine_mesh()...";
    {
      CGAL::Delaunay_mesher_2<CDT, Criteria> mesher(cdt, Criteria(0.125, 0.1));
      mesher.set_seeds(seeds.begin(), seeds.end());
      mesher.refine_mesh();
    }
    const size_type number_of_vertices2 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices() 
              << "\n\n";

    CGAL_assertion( number_of_vertices2 == number_of_vertices1bis );

    cdt = cdt2;

    std::cerr << "Meshing the triangulation with size 0.1, with\n"
              << "a loop of mesher.try_one_step_refine_mesh()...";
    size_type step = 0;
    {
      CGAL::Delaunay_mesher_2<CDT, Criteria> mesher(cdt, Criteria(0.125, 0.1));
      mesher.set_seeds(seeds.begin(), seeds.end());
      mesher.init();
      while(mesher.try_one_step_refine_mesh())
        ++step;
    }
    const size_type number_of_vertices3 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of steps: " << step << "\n\n";

    CGAL_assertion( step + inititial_number_of_vertices >= number_of_vertices3 );
    CGAL_assertion( number_of_vertices3 == number_of_vertices2 );

    cdt = cdt2;

    std::cerr << "Meshing the triangulation with size 0.1, with\n"
              << "a loop of mesher.step_by_step_refine_mesh()...";
    step = 0;
    {
      CGAL::Delaunay_mesher_2<CDT, Criteria> mesher(cdt, Criteria(0.125, 0.1));
      mesher.set_seeds(seeds.begin(), seeds.end());
      mesher.init();
      while(mesher.step_by_step_refine_mesh())
        ++step;
    }
    const size_type number_of_vertices4 = cdt.number_of_vertices();
    std::cerr << " done.\nNumber of vertices: " << cdt.number_of_vertices()
              << "\nNumber of steps: " << step << "\n\n";

    CGAL_assertion( number_of_vertices4 == number_of_vertices2 );
    CGAL_assertion( number_of_vertices4 == step + inititial_number_of_vertices );
  }
};

struct K_e_i : public CGAL::Exact_predicates_inexact_constructions_kernel {};
struct K_e_e : public CGAL::Exact_predicates_exact_constructions_kernel {};

int main()
{
  std::cerr << "TESTING WITH Exact_predicates_inexact_constructions_kernel...\n";
  Tester<K_e_i>();
  std::cerr << "TESTING WITH Exact_predicates_exact_constructions_kernel...\n";
  Tester<K_e_e>();
};
