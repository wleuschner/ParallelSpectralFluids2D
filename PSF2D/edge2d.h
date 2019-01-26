#ifndef EDGE2D_H
#define EDGE2D_H
#include "gridenums.h"

class Edge2D
{
public:
    Edge2D();
    Edge2D(unsigned int id,unsigned int v1,unsigned int v2,GridState inside);

    union
    {
        struct
        {
            unsigned int v1;
            unsigned int v2;
        };
        unsigned int v[2];
    };
    GridState inside;

    unsigned int id;
};

#endif // EDGE2D_H
