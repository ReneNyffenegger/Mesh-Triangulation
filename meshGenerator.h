#ifndef _MESHGENERATOR_H_
#define _MESHGENERATOR_H_

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <set>
#include <omp.h>
#include <simple_svg.hpp>
#include "kdtree.hpp"


using namespace svg;

class node;
class meshGenerator;

class node
{
    public:

        double pos[3];

        int index;

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
            //printf("node created: %lf %lf %lf\n", pos[0], pos[1], pos[2]);
        };


        double distanceSquared(node *p)
        {
            return( (pos[0]-p->pos[0])*(pos[0]-p->pos[0]) + (pos[1]-p->pos[1])*(pos[1]-p->pos[1]) + (pos[2]-p->pos[2])*(pos[2]-p->pos[2]));
        };

};

class triangle
{
    public:
        std::vector<node*> * nodeVector;  // pointer to the node vector
        int nodes[3];           // the indices of the 3 nodes used to construct the circle
        double pos[2];            // position of circumcenter
        double radius;            // radius of circle

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

            qsort(nodes,3, sizeof(int), intcomp);
           // printf("Constructing triangle %d %d %d\n", nodes[0],nodes[1],nodes[2]);

            constructTriangle( listOfNodes[ nodes[0]],listOfNodes[nodes[1]],listOfNodes[nodes[2]] );
        };

        void constructTriangle( node* An, node *Bn, node *Cn)
        {
            getCircle( An, Bn, Cn);
        }

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
 //           printf("Triangle created: %lf %lf, %lf (%d, %d, %d)\n", pos[0], pos[1], radius, nodes[0], nodes[1], nodes[2]);

        };

        double distanceSquared(node *A)
        {

            return( (pos[0] - A->pos[0])*(pos[0] - A->pos[0]) + (pos[1] - A->pos[1])*(pos[1] - A->pos[1]) );
        }


        static int intcomp(const void *a, const void *b)
        {
            return( *(int*)a - *(int*)b );
        }

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
        KTree<node,2> nodeTree;
        KTree<triangle,2> triangleTree;

        std::vector<node*> nodeVector;
        std::set<triangle*, triangle> triangleSet;

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

            for(int i=0;i<nodeVector.size();i++)
            {
                delete nodeVector[i];
            }


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
            if( x > 1.0 || x < 0.0 || y > 1.0 || y < 0.0) return false;

            node * n = new node(x,y,0.0);
            nodeVector.push_back(n);
            n->index = nodeVector.size()-1;
            nodeTree.insert(n,n->pos);

            return true;
        }

        void generateMesh(double maxLength)
        {
            std::vector<node*> nn;
            double pos[2];


            for(int i=0;i<nodeVector.size();i++)
            {
                nodeTree.find( nn, maxLength, nodeVector[i]->pos);

                //printf("%d, Neighbours found: %d\n", i, nn.size());
                for(int j=0 ;j < nn.size(); j++)
                {
                    for(int k=0 ;k < nn.size(); k++)
                    {
                        if( nodeVector[i] != nn[j] && nodeVector[i] != nn[k] && nn[j] != nn[k] )
                        {
                            triangle * T = new triangle( nodeVector, i, nn[j]->index, nn[k]->index);
                            int count = 0;

                            for(int l=0;l<nn.size();l++)
                            {
                                if( nn[l] != nodeVector[i] && nn[l] != nn[j] && nn[l] != nn[k] )
                                {

                                    if( T->distanceSquared(nn[l]) < T->radius || T->radius > 0.25*maxLength*maxLength){
                                        count++;
                                        break;
                                    }
                                }
                            }

                            if(count==0)
                            {
                                triangleSet.insert( T );
                            } else
                            {
                                delete T;
                            }

                        }
                    }
                }
            }

            printf("Mesh Generated\n");
            printf("    Nodes:     %d\n",(int)nodeVector.size());
            printf("    Triangles: %d\n",(int)triangleSet.size());


        }

        void print(char *outputName=NULL)
        {
            char tri[50];
            char node[50];

            if( outputName == NULL)
            {
                sprintf(tri,"mesh_triangles.txt");
                sprintf(node,"mesh_nodes.txt");
            } else {
                sprintf(tri,"%s_triangles.txt", outputName);
                sprintf(node,"%s_nodes.txt", outputName);
            }

            FILE * f = fopen(node, "w");
            for(int i=0;i<nodeVector.size();i++)
            {
                fprintf(f,"%f %f %f\n", nodeVector[i]->pos[0], nodeVector[i]->pos[1], nodeVector[i]->pos[2]);
            }
            fclose(f);

            f = fopen(tri, "w");

            std::set<triangle*>::iterator it;

            triangle *tria;
            for(it = triangleSet.begin(); it != triangleSet.end(); it++)
            {
                tria = *it;
                fprintf(f, "%d %d %d\n", tria->nodes[0], tria->nodes[1], tria->nodes[2]);
            }
            fclose(f);

            printf("Vertices printed to %s\n", node);
            printf("Triangles printed to %s\n", tri);
        }

        void drawMesh(char *outputName=NULL, double scale=1000.0)
        {

            std::set<triangle*>::iterator it;
            char fileName[50];

            if( outputName == NULL)
            {
                sprintf(fileName,"mesh.svg");
            } else {
                sprintf(fileName,"%s.svg",outputName);
            }


            Dimensions dimensions(scale, scale);
            Document doc(fileName, Layout(dimensions, Layout::BottomLeft));


            double sc = scale;
            int counter=0;

            triangle *tri;
            for(it = triangleSet.begin(); it != triangleSet.end(); it++)
            {
                //printf("drawing triangle %d\n", counter);
                counter++;

                tri = *it;

                Polyline triLine(Stroke(1, Color::Blue));



                triLine  << Point( nodeVector[tri->nodes[0]]->pos[0]*sc, nodeVector[tri->nodes[0]]->pos[1]*sc)
                         << Point( nodeVector[tri->nodes[1]]->pos[0]*sc, nodeVector[tri->nodes[1]]->pos[1]*sc)
                         << Point( nodeVector[tri->nodes[2]]->pos[0]*sc, nodeVector[tri->nodes[2]]->pos[1]*sc)
                         << Point( nodeVector[tri->nodes[0]]->pos[0]*sc, nodeVector[tri->nodes[0]]->pos[1]*sc);


//                         << Point( tri->b->x*sc, tri->b->y*sc)
                         //<< Point( tri->c->x*sc, tri->c->y*sc)
                         //<< Point( tri->a->x*sc, tri->a->y*sc);

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


