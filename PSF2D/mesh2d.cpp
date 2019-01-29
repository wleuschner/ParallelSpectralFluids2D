#include "mesh2d.h"
#include <iostream>
#include <eigen3/unsupported/Eigen/ArpackSupport>
#include <glm/gtx/rotate_vector.hpp>
#include "Spectra/GenEigsSolver.h"
#include "Spectra/MatOp/SparseGenRealShiftSolve.h"
#include "Spectra/GenEigsRealShiftSolver.h"
#include "Spectra/Util/SelectionRule.h"
#include "dec.h"
typedef Eigen::SparseMatrix<double> SparseMat;
typedef Eigen::SimplicialLDLT<SparseMat> SparseChol;
typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver <SparseMat, SparseChol> Arpack;

Mesh2D::Mesh2D()
{
}

Mesh2D::Mesh2D(unsigned int w,unsigned int h,unsigned int resolution,unsigned char* data) : Mesh2D()
{
    this->width = w;
    this->height = h;
    this->resolution = resolution;
    this->pixels = data;
}

Mesh2D::~Mesh2D()
{
}

void Mesh2D::setResolution(unsigned int res)
{
    resolution = res;
}

unsigned int Mesh2D::getWidth()
{
    return width;
}

unsigned int Mesh2D::getHeight()
{
    return height;
}

unsigned int Mesh2D::getResolution()
{
    return resolution;
}

DECMesh2D Mesh2D::voxelize()
{
    DECMesh2D decMesh(width/resolution+2,resolution);
    unsigned int iv1,iv2,iv3,iv4;

    for(unsigned int y=0;y<height/resolution+2;y++)
    {
        for(unsigned int x=0;x<width/resolution+2;x++)
        {
            if(x==0||y==0||x==width/resolution+1||y==height/resolution+1)
            {
                iv1 = y*(width/resolution+1+2)+x; //DOWN LEFT Vert
                iv2 = y*(width/resolution+1+2)+x+1; //DOWN RIGHT Vert
                iv3 = (y+1)*(width/resolution+1+2)+x; //UP LEFT vert
                iv4 = (y+1)*(width/resolution+1+2)+x+1; //UP RIGHT vert*/

                Vertex2D v1,v2,v3,v4;
                v1.pos.x = x*resolution;
                v1.pos.y = y*resolution;
                v2.pos.x = (x+1)*resolution;
                v2.pos.y = y*resolution;
                v3.pos.x = x*resolution;
                v3.pos.y = (y+1)*resolution;
                v4.pos.x = (x+1)*resolution;
                v4.pos.y = (y+1)*resolution;

                vertex[iv1] = v1;
                vertex[iv2] = v2;
                vertex[iv3] = v3;
                vertex[iv4] = v4;

                decMesh.addFace(Face2D(y*(width/resolution+2)+x,iv1,iv3,iv4,iv2,GridState::OUTSIDE));
            }
            else
            {
                unsigned int ix=x-1;
                unsigned int iy=y-1;
                unsigned char* start = &pixels[iy*width*resolution+ix*resolution];
                iv1 = y*(width/resolution+1+2)+x; //DOWN LEFT Vert
                iv2 = y*(width/resolution+1+2)+x+1; //DOWN RIGHT Vert
                iv3 = (y+1)*(width/resolution+1+2)+x; //UP LEFT vert
                iv4 = (y+1)*(width/resolution+1+2)+x+1; //UP RIGHT vert*/

                bool inside = checkVoxel(start);
                Vertex2D v1,v2,v3,v4;
                v1.pos.x = x*resolution;
                v1.pos.y = y*resolution;
                v2.pos.x = (x+1)*resolution;
                v2.pos.y = y*resolution;
                v3.pos.x = x*resolution;
                v3.pos.y = (y+1)*resolution;
                v4.pos.x = (x+1)*resolution;
                v4.pos.y = (y+1)*resolution;

                vertex[iv1] = v1;
                vertex[iv2] = v2;
                vertex[iv3] = v3;
                vertex[iv4] = v4;

                if(inside)
                {
                    decMesh.addFace(Face2D(y*(width/resolution+2)+x,iv1,iv3,iv4,iv2,GridState::INSIDE));
                }
                else
                {
                    decMesh.addFace(Face2D(y*(width/resolution+2)+x,iv1,iv3,iv4,iv2,GridState::OUTSIDE));
                }
            }
        }
    }
    return decMesh;
}

bool Mesh2D::checkVoxel(unsigned char* offs)
{
    for(unsigned int ry=0;ry<resolution;ry++)
    {
        for(unsigned int rx=0;rx<resolution;rx++)
        {
            if(offs[ry*width+rx]==0)
            {
                return true;
            }
        }
    }
    return false;
}
