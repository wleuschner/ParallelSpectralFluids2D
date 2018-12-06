#include "edge2d.h"
#include "mesh2d.h"
#include <cmath>

Edge2D::Edge2D()
{
    v1 = -1;
    v2 = -1;
}

Edge2D::Edge2D(unsigned int v1,unsigned int v2) : Edge2D()
{
    this->v1 = v1;
    this->v2 = v2;
}
