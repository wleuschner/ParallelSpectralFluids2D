#ifndef MESH2D_H
#define MESH2D_H
#include "face2d.h"
#include "edge2d.h"
#include "vertex2d.h"
#include "decmesh2d.h"
#include <map>
#include <set>
#include <Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

class Mesh2D
{
public:
    Mesh2D();
    Mesh2D(unsigned int w,unsigned int h,unsigned int resolution,unsigned char* data);
    ~Mesh2D();

    void setResolution(unsigned int res);

    unsigned int getWidth();
    unsigned int getHeight();
    unsigned int getResolution();
    DECMesh2D voxelize();

private:
    bool checkVoxel(unsigned char* offs);

    unsigned int width;
    unsigned int height;
    unsigned int resolution;
    unsigned char* pixels;

public:
    std::map<unsigned int,Vertex2D> vertex;
};

#endif // MESH2D_H
