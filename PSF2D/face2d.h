#ifndef FACE2D_H
#define FACE2D_H
#include"edge2d.h"

class Mesh2D;
class Edge2D;

class Face2D
{
public:
    Face2D();
    Face2D(int v1,int v2,int v3,int v4);

    union
    {
        struct
        {
            int v1;
            int v2;
            int v3;
            int v4;
        };
        int v[4];
    };
};

#endif // FACE2D_H
