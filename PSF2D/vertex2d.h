#ifndef VERTEX2D_H
#define VERTEX2D_H
#include "edge2d.h"
#include <glm/glm.hpp>

class Mesh2D;
class Edge2D;

class Vertex2D
{
public:
    Vertex2D();
    Vertex2D(unsigned int x,unsigned int y);

    glm::vec2 pos;
};

#endif // VERTEX2D_H
