#ifndef EDGE2D_H
#define EDGE2D_H
#include "vertex2d.h"
#include "face2d.h"

class Mesh2D;
class Face2D;

class Edge2D
{
public:
    Edge2D();
    Edge2D(unsigned int v1,unsigned int v2);
    float area;

    union
    {
        struct
        {
            unsigned int v1;
            unsigned int v2;
        };
        unsigned int v[2];
    };

    int faces[2];
};

#endif // EDGE2D_H
