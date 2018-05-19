#ifndef BENEATH_BEYOND_3
#define BENEATH_BEYOND_3

#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Timer.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

using namespace std;

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
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
        typedef typename Facet::Plane_3 Plane;
        return Plane(h->vertex()->point(),
                     h->next()->vertex()->point(),
                     h->next()->next()->vertex()->point());
    }
};


template <typename Iterator, typename objType>
void beneath_beyond_3(Iterator begin, Iterator end, objType &obj) {
    // sort points lexicographical
    sort(begin, end);
    
    Iterator itr = begin;
    Iterator initIters[4];
    
    for (int i = 0; i < 4; ++i) {
        initIters[i] = itr;
        ++itr;
    }
    
    // create initial tetrahedron
    obj.make_tetrahedron(*initIters[0], *initIters[1], *initIters[2], *initIters[3]);
     
    for (; itr != end; ++itr) {
        for (FacetIterator facetItr = obj.facets_begin(); facetItr != obj.facets_end(); ++facetItr) {
            HalfEdgeHandle edge = facetItr->halfedge();
            Point3 pnt1 = edge->vertex()->point();
            Point3 pnt2 = edge->next()->vertex()->point();
            Point3 pnt3 = edge->next()->next()->vertex()->point();
            
            //cout << pnt1 << " " << pnt2 << " " << pnt3 << " compare with " << *itr << " == ";
            //cout << CGAL::orientation(pnt1, pnt2, pnt3, *itr) << endl; 

            if (CGAL::orientation(pnt1, pnt2, pnt3, *itr) == 1) {
                obj.erase_facet(edge);
                make_triangle();
            }
        }
    }
}

#endif
