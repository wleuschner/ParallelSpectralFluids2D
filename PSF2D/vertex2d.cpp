#include "vertex2d.h"

Vertex2D::Vertex2D()
{
    inside = GridState::UNINITIALIZED;
}

Vertex2D::Vertex2D(unsigned int v,GridState inside) : Vertex2D()
{
    this->inside = inside;
    this->id = v;
}
