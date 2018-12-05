#include "face2d.h"
#include "mesh2d.h"

Face2D::Face2D()
{
    area = 0.0f;
    e1 = -1;
    e2 = -1;
    e3 = -1;
    e4 = -1;
}

Face2D::Face2D(int e1,int e2,int e3,int e4) : Face2D()
{
    this->e1 = e1;
    this->e2 = e2;
    this->e3 = e3;
    this->e4 = e4;
}
