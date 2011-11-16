#ifndef _MESHGENERATOR_H_
#define _MESHGENERATOR_H_

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <set>
#include <omp.h>
#include <string>
#include "simple_svg.hpp"
#include "kdtree.hpp"


using namespace svg;

class node;
class meshGenerator;


/*==========================================================
A node class in 3 dimensions (only 2 are used so far)
==========================================================*/
class node
{
    public:

        double pos[3];                          // The cartesian coordinates of the node
        int index;                              // The index of the node.  This is automatically modified. Do not play with thiss

        node& operator=(const node & n)
        {
            pos[0] = n.pos[0];
            pos[1] = n.pos[1];
            pos[2] = n.pos[2];
        }

        node(double X=0.0, double Y=0.0, double Z=0.0)
        {
            pos[0] = X;
            pos[1] = Y;
            pos[2] = Z;
        };

        // Calculates the square of the distance between this node and another node
        double distanceSquared(node *p)
        {
            return( (pos[0]-p->pos[0])*(pos[0]-p->pos[0]) + (pos[1]-p->pos[1])*(pos[1]-p->pos[1]) + (pos[2]-p->pos[2])*(pos[2]-p->pos[2]));
        };

};


/*==========================================================
Triangle class

Information about a triangle in the mesh.
==========================================================*/
class triangle
{
    public:
        std::vector<node*> * nodeVector;                // pointer to the node vector
        int nodes[3];                                   // the indices of the 3 nodes used to construct the circle
        double pos[2];                                  // position of circumcenter
        double radius;                                  // radius of circle

        triangle()
        {
                nodeVector = NULL;
        }

        triangle& operator=(const triangle & n)
        {
            nodeVector = n.nodeVector;
            nodes[0] = n.nodes[0];
            nodes[1] = n.nodes[1];
            nodes[2] = n.nodes[2];
            pos[0] = n.pos[0];
            pos[1] = n.pos[1];
            radius = n.radius;
        }

        triangle(std::vector<node*> &listOfNodes, int i, int j, int k)
        {
            nodes[0] = i;
            nodes[1] = j;
            nodes[2] = k;

            nodeVector = &listOfNodes;

            // sort the node indices so they are from smallest index to largest index
            qsort(nodes,3, sizeof(int), intcomp);

            constructTriangle( listOfNodes[ nodes[0]],listOfNodes[nodes[1]],listOfNodes[nodes[2]] );
        };

        void constructTriangle( node* An, node *Bn, node *Cn)
        {
            getCircle( An, Bn, Cn);
        }

        /*==============================================================
        Given the 3 nodes of the triangle, calculate the circumcenter
        and the radius of the circle
        ==============================================================*/
        void getCircle( node* A, node *B, node *C )
        {

            double Bx = B->pos[0] - A->pos[0];
            double Cx = C->pos[0] - A->pos[0];

            double By = B->pos[1] - A->pos[1];
            double Cy = C->pos[1] - A->pos[1];

            double D = 2*(Bx*Cy - By*Cx);

            double nx = ( Cy*(Bx*Bx + By*By) - By*(Cx*Cx + Cy*Cy) ) / D;
            double ny = ( Bx*(Cx*Cx + Cy*Cy) - Cx*(Bx*Bx + By*By) ) / D;

            pos[0] = nx+A->pos[0];
            pos[1] = ny+A->pos[1];


            radius = (pos[0] - A->pos[0])*(pos[0] - A->pos[0]) + (pos[1] - A->pos[1])*(pos[1] - A->pos[1]);

        };

        /*==============================================================
        Square of the distance between a node and the circumcenter of
        the circle created by the triangle
        ==============================================================*/
        double distanceSquared(node *A)
        {

            return( (pos[0] - A->pos[0])*(pos[0] - A->pos[0]) + (pos[1] - A->pos[1])*(pos[1] - A->pos[1]) );
        }

        // Compares integers for the qsort algorithm
        static int intcomp(const void *a, const void *b)
        {
            return( *(int*)a - *(int*)b );
        }

        // compairson operator for the std::vector<triangle*> function
        bool operator()(triangle* a, triangle *b)
        {
            if( a->nodes[0] < b->nodes[0]) return(true);
            if( a->nodes[1] < b->nodes[1]) return(true);
            if( a->nodes[2] < b->nodes[2]) return(true);

            return(false);
        }
};



class meshGenerator
{
    private:
        KTree<node,2>                   nodeTree;                       // A KD-Tree of nodes. This is used to find the neasrest neighbours
        KTree<triangle,2>               triangleTree;                   // A KD-Tree for the triangles, This isn't used it. Maybe in the future.

        std::vector<node*>              nodeVector;                     // A vector of all the nodes in the mesh
        std::set<triangle*, triangle>   triangleSet;                    // A set of all the triangles. Cannot contain duplicate triangles

        double maxX;
        double minX;

        double maxY;
        double minY;

    public:

        meshGenerator()
        {
            maxX = maxY = -999999.9;
            minX = minY = 999999.9;
        };

        ~meshGenerator()
        {

            // Delete all the nodes that have been created
            for(int i=0;i<nodeVector.size();i++)
            {
                delete nodeVector[i];
            }


            // delete all the triangles that have been created
            std::set<triangle*>::iterator it;
            triangle *tria;
            for(it = triangleSet.begin(); it != triangleSet.end(); it++)
            {
                tria = *it;
                delete tria;
            }



        }

        bool insertNode(double x, double y)
        {
            // Make sure the node we are inserting is in the range (0,0 <= (x,y) < (1,1)
            if( x > 1.0 || x < 0.0 || y > 1.0 || y < 0.0) return false;

            node * n = new node(x,y,0.0);

            // add it to the node vector
            nodeVector.push_back(n);

            // set the index of the new node
            n->index = nodeVector.size()-1;

            // insert the node into the KD-tree
            nodeTree.insert(n,n->pos);

            return true;
        }


        /*=============================================================================
        Generates the mesh out of the nodes that are already in the vector.

        Outline of the algorithm
           -For every 3 points, construct a circle that passes through all three points
           -If there are no other points that lie within that circle, then construct
            a triangle with those three points.

        Input paramter: maxLength - the maximum length of the side of a triangle
                                  - this is not strictly enforced.
                                  - This parameter is to increase the speed of the
                                    algorithm.
        *============================================================================*/
        void generateMesh(double maxLength)
        {
            std::vector<node*>              nn;                // a vector of nearest neighbours

            // For all nodes in the node vector
            for(int i=0;i<nodeVector.size();i++)
            {
                // find all the nearest neighbours that are within maxLength and store them in the vector nn.
                nodeTree.find( nn, maxLength, nodeVector[i]->pos);

                // loop through all the pairs in the nearest neighbours vector
                for(int j=0 ;j < nn.size(); j++)
                for(int k=0 ;k < nn.size(); k++)
                {
                    // make sure that all three nodes are different.
                    if( nodeVector[i] != nn[j] && nodeVector[i] != nn[k] && nn[j] != nn[k] )
                    {
                        //construct the triangle with those three nodes
                        triangle * T = new triangle( nodeVector, i, nn[j]->index, nn[k]->index);
                        int count = 0;

                        //loop through all the nearest neighbours,
                        for(int l=0;l<nn.size();l++)
                        {
                            // making sure we do not include the vertices of the triangle
                            if( nn[l] != nodeVector[i] && nn[l] != nn[j] && nn[l] != nn[k] )
                            {
                                // and check if the node is within the circle of the triangle.
                                if( T->distanceSquared(nn[l]) < T->radius || T->radius > 0.25*maxLength*maxLength){
                                    count++;
                                    break;
                                }
                            }
                        }

                        // if no other nodes are within the circle
                        if(count==0)
                        {
                            // add the triangle to the list
                            triangleSet.insert( T );
                        } else
                        {
                            //otherwise, delete the triangle
                            delete T;
                        }

                    }
                }

            }

            printf("Mesh Generated\n");
            printf("    Nodes:     %d\n",(int)nodeVector.size());
            printf("    Triangles: %d\n",(int)triangleSet.size());


        }

        void print(std::string outputName)
        {
	    std::string tri;
	    std::string node;

	    tri = outputName + "_triangles.txt";
	    node = outputName + "_nodes.txt";

            FILE * f = fopen(node.c_str(), "w");
            for(int i=0;i<nodeVector.size();i++)
            {
                fprintf(f,"%f %f %f\n", nodeVector[i]->pos[0], nodeVector[i]->pos[1], nodeVector[i]->pos[2]);
            }
            fclose(f);

            f = fopen(tri.c_str(), "w");

            std::set<triangle*>::iterator it;

            triangle *tria;
            for(it = triangleSet.begin(); it != triangleSet.end(); it++)
            {
                tria = *it;
                fprintf(f, "%d %d %d\n", tria->nodes[0], tria->nodes[1], tria->nodes[2]);
            }
            fclose(f);

            printf("Vertices printed to %s\n", node.c_str());
            printf("Triangles printed to %s\n", tri.c_str());
        }

        void drawMesh(std::string outputName, double scale=1000.0)
        {

            std::set<triangle*>::iterator it;
            char fileName[50];

            if( outputName == "")
            {
                sprintf(fileName,"mesh.svg");
            } else {
                sprintf(fileName,"%s.svg",outputName.c_str());
            }


            Dimensions dimensions(scale, scale);
            Document doc(fileName, Layout(dimensions, Layout::BottomLeft));


            double sc = scale;
            int counter=0;

            triangle *tri;
            for(it = triangleSet.begin(); it != triangleSet.end(); it++)
            {

                counter++;

                tri = *it;

                Polyline triLine(Stroke(1, Color::Blue));



                triLine  << Point( nodeVector[tri->nodes[0]]->pos[0]*sc, nodeVector[tri->nodes[0]]->pos[1]*sc)
                         << Point( nodeVector[tri->nodes[1]]->pos[0]*sc, nodeVector[tri->nodes[1]]->pos[1]*sc)
                         << Point( nodeVector[tri->nodes[2]]->pos[0]*sc, nodeVector[tri->nodes[2]]->pos[1]*sc)
                         << Point( nodeVector[tri->nodes[0]]->pos[0]*sc, nodeVector[tri->nodes[0]]->pos[1]*sc);


                doc << triLine;

            }

            for(int i=0;i<nodeVector.size();i++)
            {
                doc << Circle(Point( nodeVector[i]->pos[0]*sc, nodeVector[i]->pos[1]*sc), 10, Fill(Color(100, 200, 120)), Stroke(0, Color(200, 250, 150)));
            }

            doc.save();
            printf("Mesh printed to %s\n", fileName);

        }
};

/*
class meshGenerator1
{
    private:
        void * kd;
        void * kdTriangle;

            KTree<node,2> nodeTree;
        KTree<triangle,2> triangleTree;

        std::vector<triangle*> triangleList;

        std::vector<node*> listOfNodes;

        //std::vector<int*> triangleList;

        std::set<triangle*, triangleCompare> triangleSet;

    public:
        meshGenerator()
        {
            createKDTree();
        }

        void test()
        {

            node center;
            #define co (3.14159/180.0)
            node a( cos(0)     , sin(0)      );
            node b( cos(180*co), sin(180*co) );
            node c( cos(90*co) , sin(90*co)  );

            center = getCenterOfCircle( &a, &b, &c);

            printf("Center: %f %f,  Radius = %f\n", center.x, center.y, center.distanceSquared(&a));

        }
        ~meshGenerator()
        {
            for(int i=0;i<listOfNodes.size();i++)
            {
                delete listOfNodes[i];
            }
            for(int i=0;i<triangleList.size();i++)
            {
                delete []  triangleList[i];
            }
            kd_free( (kdtree*)kd);
        }

        void createKDTree()
        {
            kd = kd_create( 3 );
            kdTriangle = kd_create( 3 );
        }

        void insertNode( double x, double y, double z=0.0)
        {
            node * n = new node( x, y, z);
            insertNode(n);
        }


        void generateMesh(double largestEdgeSize)
        {
            std::vector<node*> nn;


            for(int i=0;i<listOfNodes.size();i++)
            {
                nodeTree.find(largestEdgeSize, );

               // #pragma omp parallel for
                for(int k=0;k< nn.size();k++)
                {
                    for(int l=0; l<nn.size();l++)
                    {
                        //printf("Generating new triangle\n");


                        if( listOfNodes[i] != nn[k] && listOfNodes[i] != nn[l] && nn[k] != nn[l] )
                        {
                            node center;
                            center = getCenterOfCircle( listOfNodes[i], nn[k], nn[l]);

                            double r2 = center.distanceSquared( listOfNodes[i] );

                            if( r2 > 0.25*largestEdgeSize*largestEdgeSize)
                            {
                                    continue;
                                   // printf("found large triangle\n");
                            }

                            int inCenter = 0;

                            for(int m=0; m < nn.size(); m++)
                            {
                                if( nn[m] != listOfNodes[i] && nn[m] != nn[l] && nn[m] != nn[k])
                                {
                                    double d2 = center.distanceSquared( nn[m] );

                                    if( d2 < r2)
                                    {
                                        inCenter++;
                                        m=nn.size();
                                        break;
                                    }
                                }
                            }

                            if( inCenter == 0)
                            {

                                insertTriangle( new triangle( listOfNodes, listOfNodes[i]->index, nn[k]->index,  nn[l]->index ) );

                                //triangleList.push_back(nod);
                                //triangleSet.insert(nod);

                            } else {
                               // printf("Triangle rejected\n");
                            }

                        }
                    }
                }
                //printf("NODE COMPLETE\n");
            }
        }

        static node getCenterOfCircle( node* A, node *B, node *C )
        {

            double Bx = B->x - A->x;
            double Cx = C->x - A->x;
            double By = B->y - A->y;
            double Cy = C->y - A->y;

            double D = 2*(Bx*Cy - By*Cx);

            double nx = ( Cy*(Bx*Bx + By*By) - By*(Cx*Cx + Cy*Cy) ) / D;
            double ny = ( Bx*(Cx*Cx + Cy*Cy) - Cx*(Bx*Bx + By*By) ) / D;

            node n(nx+A->x,ny+A->y, 0.0 );

            return(n);
        };

        void drawMesh()
        {

            std::set<int*>::iterator it;
            int *tri;

            double sc = 1000;
            Dimensions dimensions(1000, 1000);
            Document doc("my_svg2.svg", Layout(dimensions, Layout::BottomLeft));



            int counter=0;

            for(it = triangleSet.begin(); it != triangleSet.end(); it++)
            {
                printf("drawing triangle %d\n", counter);
                counter++;

                tri = *it;
                int a = tri[0];
                int b = tri[1];
                int c = tri[2];

                Polyline triLine(Stroke(1, Color::Blue));
                triLine  << Point( listOfNodes[a]->x*sc, listOfNodes[a]->y*sc)
                                    << Point( listOfNodes[b]->x*sc, listOfNodes[b]->y*sc)
                                    << Point( listOfNodes[c]->x*sc, listOfNodes[c]->y*sc)
                                    << Point( listOfNodes[a]->x*sc, listOfNodes[a]->y*sc);
                doc << triLine;

            }

            for(int i=0;i<listOfNodes.size();i++)
            {
                doc << Circle(Point( listOfNodes[i]->x*sc, listOfNodes[i]->y*sc), 10, Fill(Color(100, 200, 120)), Stroke(0, Color(200, 250, 150)));
            }

            doc.save();


        }

        void octavePrint()
        {
            FILE * f = fopen("output.m", "w");
            std::set<int*>::iterator it;



            fprintf(f,"nodes = [");
            for(int i=0;i<listOfNodes.size()-1;i++)
            {
                fprintf(f,"%f,%f;", listOfNodes[i]->x, listOfNodes[i]->y);
            }
            fprintf(f,"%f,%f];", listOfNodes[ listOfNodes.size()-1]->x, listOfNodes[listOfNodes.size()-1]->y);




            int *tri;
            fprintf(f,"\n\ntri = [");

            for(it = triangleSet.begin(); it != triangleSet.end(); it++)
            {
                tri = *it;
                printf("%d,%d,%d\n", tri[0], tri[1], tri[2]);
                fprintf(f,"%d,%d,%d\n", tri[0], tri[1], tri[2]);

            }
            fprintf(f,"];");

            printf("\n\nPrinted to file\n");
            printf("   Nodes:     %d\n", (int)listOfNodes.size());
            printf("   Triangles: %d\n", (int)triangleSet.size());

            fclose(f);
        }

        private:

            static int intcomp(const void *a, const void *b)
            {
                return( *(int*)a - *(int*)b );
            }


            void insertTriangle(triangle * n)
            {
                double pos[2] = {n->x, n->y};

                if( !nodeTree.insert( n, pos) )
                {
                    listOfNodes.push_back(n);
                    n->index = listOfNodes.size()-1;
                }
            }

            void insertNode(node * n)
            {
                double pos[2] = {n->x, n->y};

                if( !nodeTree.insert( n, pos) )
                {
                    listOfNodes.push_back(n);
                    n->index = listOfNodes.size()-1;
                }

            }



};

*/

#endif // _MESHGENERATOR_H_


