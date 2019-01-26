#include "edge2d.h"
#include <cmath>

Edge2D::Edge2D()
{
    inside = GridState::UNINITIALIZED;
}

Edge2D::Edge2D(unsigned int id,unsigned int v1,unsigned int v2,GridState inside) : Edge2D()
{
    this->inside = inside;
    this->id = id;
    this->v1 = v1;
    this->v2 = v2;
}
