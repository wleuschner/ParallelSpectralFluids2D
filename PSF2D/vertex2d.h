#ifndef VERTEX2D_H
#define VERTEX2D_H
#include <glm/glm.hpp>
#include "gridenums.h"

class Vertex2D
{
public:
    Vertex2D();
    Vertex2D(unsigned int v,GridState inside);

    union
    {
        struct
        {
            int e1;
            int e2;
            int e3;
            int e4;
        };
        int e[4];
    };

    GridState inside;
    unsigned int id;

    glm::dvec2 pos;
};

#endif // VERTEX2D_H
