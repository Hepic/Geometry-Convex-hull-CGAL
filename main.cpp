#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "beneath_beyond_3.h"

using namespace std;

void printInformation(Polyhedron3 polyHed, CGAL::Timer timer) {
    // print vertices
    for (VertexIterator itr = polyHed.vertices_begin(); itr != polyHed.vertices_end(); ++itr) {
        cout << itr->point() << endl;
    }
    
    // print facets as plane equations
    transform(polyHed.facets_begin(), polyHed.facets_end(), polyHed.planes_begin(), planeEquation());
    copy(polyHed.planes_begin(), polyHed.planes_end(), ostream_iterator<Plane3>(cout, "\n"));
    
    // print time of calculation of convex hull
    cout << "Time passed: " << timer.time() << " seconds" << endl;
}


int main(int argc, char *argv[]) {
    CGAL::set_pretty_mode(std::cout);
    
    int N;
    ifstream inputFile;
    vector<Point3> points;
    Polyhedron3 polyHed;
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
    
    // calculate convex hull with algorithm from question 1
    timer.start();

    for (int i = 0; i < N; ++i) {
        convHull.insert(points[i]);
    }
    
    CGAL::convex_hull_d_to_polyhedron_3(convHull, polyHed);
    timer.stop();

    cout << "Convex hull from question 1\n---------------------------\n";
    printInformation(polyHed, timer);
    
    // calculate convex hull with algorithm from question 2
    timer.start();
    polyHed.clear();
    
    // returns true only if there is no degeneracy
    if (beneath_beyond_3(points.begin(), points.end(), polyHed)) {
        timer.stop();
        
        cout << "\nConvex hull from question 2\n---------------------------\n";
        printInformation(polyHed, timer);
    }
    
    return 0;
}
