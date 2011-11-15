#include "meshGenerator.h"
#include <simple_svg.hpp>
#include <ting/Socket.hpp>
#include <time.h>

#define N 250
#define RANDOM ( (double)rand() / (double) RAND_MAX )

int main(int argc, char **argv)
{
    meshGenerator M;
    srand( time(NULL) );

    // Generate some random points
    //
    // Points must be in the range:
    // (0.0) < (x,y) < (1.0,1.0)
    for(int i=0;i<N;i++)
    {
        double x = RANDOM;
        double y = RANDOM;
        double th = RANDOM*2*3.14159;
        M.insertNode( x,y);
    }

    // Generate a mesh making sure that the max length of
    // a triangle's edige is less than 0.25;
    M.generateMesh(0.25);

    // draw the mesh to an svg file that has dimensions 1000 x 1000
    M.drawMesh("mesh", 1000.0);

    // print the mesh to a text file (mesh_triangles.txt) and (mesh_vertices.txt)
    M.print("mesh");

    return 0;
}
