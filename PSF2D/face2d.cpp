#include "face2d.h"
#include "mesh2d.h"

Face2D::Face2D()
{
    v1 = -1;
    v2 = -1;
    v3 = -1;
    v4 = -1;
}

Face2D::Face2D(int e1,int e2,int e3,int e4) : Face2D()
{
    this->v1 = e1;
    this->v2 = e2;
    this->v3 = e3;
    this->v4 = e4;
}
