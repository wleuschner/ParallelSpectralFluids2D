#ifndef FACE2D_H
#define FACE2D_H
#include "gridenums.h"

class Face2D
{
public:
    Face2D();
    Face2D(unsigned int id,unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,GridState inside);

    union
    {
        struct
        {
            unsigned int v1;
            unsigned int v2;
            unsigned int v3;
            unsigned int v4;
        };
        unsigned int v[4];
    };

    union
    {
        struct
        {
            unsigned int e1;
            unsigned int e2;
            unsigned int e3;
            unsigned int e4;
        };
        unsigned int e[4];
    };
    GridState inside;

    unsigned int id;
};

#endif // FACE2D_H
