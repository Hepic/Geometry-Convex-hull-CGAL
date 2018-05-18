#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Timer.h>

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point3;
typedef Kernel::Plane_3 Plane3;
typedef CGAL::Creator_uniform_3<double, Point3> Creator;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron3;
typedef Polyhedron3::Vertex_iterator VertexIterator;
typedef Polyhedron3::Facet_iterator FacetIterator;
typedef Polyhedron3::Halfedge_handle HalfEdgeHandle;
typedef CGAL::Convex_hull_d_traits_3<Kernel> HullTraits3;
typedef CGAL::Convex_hull_d<HullTraits3> ConvexHull3;

struct planeEquation {
    template <class Facet>
    typename Facet::Plane_3 operator()(Facet &f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        typedef typename Facet::Plane_3  Plane;
        return Plane(h->vertex()->point(),
                     h->next()->vertex()->point(),
                     h->next()->next()->vertex()->point());
    }
};


int main(int argc, char *argv[]) {
    CGAL::set_pretty_mode(std::cout);
    
    int N;
    ifstream inputFile;
    vector<Point3> points;
    ConvexHull3 convHull(3);
    CGAL::Timer timer;

    if (argc < 2) {
        cerr << "Not enough arguments" << endl;
        return 0;
    } 
    
    if (!strcmp(argv[1], "-generate")) {
        N = atoi(argv[2]);
        CGAL::Random_points_in_sphere_3<Point3, Creator> generator(150.0);
        CGAL::copy_n(generator, N, back_inserter(points));
    } else {
        inputFile.open(argv[1]);

        if (inputFile.is_open()) {
            inputFile >> N;
            
            for (int i = 0; i < N; ++i) {
                Point3 pnt;
                
                // read points from input file
                inputFile >> pnt;
                points.push_back(pnt);
            }
        }

        inputFile.close();
    }
    
    // calculate convex hull
    timer.start();

    for (int i = 0; i < N; ++i) {
        convHull.insert(points[i]);
    }
    
    Polyhedron3 polyHed;
    CGAL::convex_hull_d_to_polyhedron_3(convHull, polyHed);

    timer.stop();

    // print vertices
    for (VertexIterator itr = polyHed.vertices_begin(); itr != polyHed.vertices_end(); ++itr) {
        cout << itr->point() << endl;
    }
    
    // print facets as plane equations
    transform(polyHed.facets_begin(), polyHed.facets_end(), polyHed.planes_begin(), planeEquation());
    copy(polyHed.planes_begin(), polyHed.planes_end(), ostream_iterator<Plane3>(cout, "\n"));
    
    // print time of calculation of convex hull
    cout << "Time passed: " << timer.time() << " seconds" << endl;

    return 0;
}
