#include "edge2d.h"
#include "mesh2d.h"
#include <cmath>

Edge2D::Edge2D()
{
    faces[0] = -1;
    faces[1] = -1;

    v1 = -1;
    v2 = -1;

    area = 0.0f;
}

Edge2D::Edge2D(unsigned int v1,unsigned int v2) : Edge2D()
{
    this->v1 = v1;
    this->v2 = v2;
}
