#include "vertex2d.h"
#include "mesh2d.h"

Vertex2D::Vertex2D()
{
    edges[0] = -1;
    edges[1] = -1;
    edges[2] = -1;
    edges[3] = -1;

    area = 0.0f;
}

Vertex2D::Vertex2D(unsigned int x,unsigned int y) : Vertex2D()
{
    this->pos.x = x;
    this->pos.y = y;
}
