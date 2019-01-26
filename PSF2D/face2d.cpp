#include "face2d.h"

Face2D::Face2D()
{
    inside = GridState::UNINITIALIZED;
}

Face2D::Face2D(unsigned int id,unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,GridState inside) : Face2D()
{
    this->inside = inside;
    this->id = id;

    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
    this->v4 = v4;
}
