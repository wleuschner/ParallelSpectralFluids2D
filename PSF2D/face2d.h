#ifndef FACE2D_H
#define FACE2D_H
#include"edge2d.h"

class Mesh2D;
class Edge2D;

class Face2D
{
public:
    Face2D();
    Face2D(int e1,int e2,int e3,int e4);

    float area;

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
};

#endif // FACE2D_H
