// Delaunay Triangulation of a set of 3D points in the xy-plane.
// (Terrain triangulation)

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>


using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_3<Kernel>  Point_3;

typedef CGAL::Projection_traits_xy_3<Kernel>               Gt;

typedef  CGAL::Triangulation_2<Gt>                 Triangulation;
typedef  CGAL::Delaunay_triangulation_2<Gt>        Delaunay_triangulation;

bool  verbose      = false;
bool  delaunay     = false;
bool  incr         = false;

// main function with standard unix commandline arguments
// ------------------------------------------------------
int main( int argc, char **argv) {
    int n = 0; // number of filenames
    char *filename[2];
    bool help = false;
    for (int i = 1; i < argc; i++) { // check commandline options
        if ( strcmp( "-v", argv[i]) == 0)
            verbose = true;
        else if ( strcmp( "-delaunay", argv[i]) == 0)
            delaunay = true;
        else if ( strcmp( "-incr", argv[i]) == 0)
            incr = true;
        else if ( (strcmp( "-h", argv[i]) == 0) ||
                  (strcmp( "-help", argv[i]) == 0))
            help = true;
        else if ( n < 2 ) {
            filename[ n++] = argv[i];
        } else {
            ++n;
            break;
        }
    }
    if ((n > 2) || help) {
        if ( ! help)
            cerr << "Error: in parameter list" << endl;
        cerr << "Usage: " << argv[0] << " [<options>] [<infile> [<outfile>]]"
             << endl;
        cerr << "       Terrain triangulation in the xy-plane." << endl;
        cerr << "       -delaunay  Delaunay triangulation (default)." << endl;
        cerr << "       -incr      Incremental insertion (no flips)." << endl;
        cerr << "       -v         verbose." << endl;
        exit( ! help);
    }

    CGAL::Verbose_ostream vout( verbose);
    vout << argv[0] << ": verbosity on." << endl;

    const char*  iname = "cin";
    istream*     p_in  = &cin;
    ifstream     in;
    if ( n > 0) {
        in.open( filename[0]);
        p_in = &in;
        iname = filename[0];
    }
    if ( !*p_in) {
        cerr << argv[0] << ": error: cannot open file '" << iname
         << "' for reading." <<endl;
        exit( 1);
    }

    CGAL::File_scanner_OFF scanner( * p_in, true);
    if ( !*p_in)
        exit( 1);

    const char*  oname = "cout";
    ostream*     p_out = &cout;
    ofstream     out;
    if ( n > 1) {
        out.open( filename[1]);
        p_out = &out;
        oname = filename[1];
    }
    if ( !*p_out) {
        cerr << argv[0] << ": error: cannot open file '"<< oname
             << "' for writing." <<endl;
        exit( 1);
    }

    if ( delaunay || ! incr) {
        Delaunay_triangulation triang;
        vout << "Scanning and triangulating ..." << endl;
        for ( std::size_t j = 0; j < scanner.size_of_vertices(); j++) {
            double x, y, z;
            scanner.scan_vertex( x, y, z);
            Point_3 p( x, y, z);
            triang.insert( p);
        }
        vout << "    .... done." << endl;
        vout << "write_triangulation( " << oname << ") ...." << endl;
        CGAL::IO::write_OFF(*p_out, triang);
        vout << "    .... done." << endl;
    } else {
        Triangulation triang;
        vout << "Scanning and triangulating ..." << endl;
        for ( std::size_t j = 0; j < scanner.size_of_vertices(); j++) {
            double x, y, z;
            scanner.scan_vertex( x, y, z);
            Point_3 p( x, y, z);
            triang.insert( p);
        }
        vout << "    .... done." << endl;
        vout << "write_triangulation( " << oname << ") ...." << endl;
        CGAL::IO::write_OFF(*p_out, triang);
        vout << "    .... done." << endl;
    }
    if ( !*p_in) {
        cerr << argv[0] << " read error: while reading file '"<< iname << "'."
             << endl;
        exit( 1);
    }
    if ( !*p_out) {
        cerr << argv[0] << " write error: while writing file '"<< oname << "'."
             << endl;
        exit( 1);
    }

    return 0;
}
