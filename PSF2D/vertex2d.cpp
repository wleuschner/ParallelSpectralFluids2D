#include "vertex2d.h"
#include "mesh2d.h"

Vertex2D::Vertex2D()
{

}

Vertex2D::Vertex2D(unsigned int x,unsigned int y) : Vertex2D()
{
    this->pos.x = x;
    this->pos.y = y;
}
