#ifndef MESH2D_H
#define MESH2D_H
#include "face2d.h"
#include "edge2d.h"
#include "vertex2d.h"
#include <map>
#include <set>
#include <Eigen/Eigen>

class Mesh2D
{
public:
    Mesh2D();
    Mesh2D(unsigned int w,unsigned int h,unsigned int resX,unsigned int resY,unsigned char* data);
    ~Mesh2D();
    unsigned int getWidth();
    unsigned int getHeight();
    bool* getVoxelMap();
    unsigned int getResX();
    unsigned int getResY();
    unsigned int getNumVertices();
    unsigned int getNumEdges();
    unsigned int getNumFaces();
    Eigen::VectorXf getVelocityField();
    Eigen::SparseMatrix<float> buildLaplace();

    friend Eigen::SparseMatrix<float> hodge2(Mesh2D& mesh,bool inverse);
    friend Eigen::SparseMatrix<float> hodge1(Mesh2D& mesh,bool inverse);
    friend Eigen::SparseMatrix<float> hodge0(Mesh2D& mesh,bool inverse);

    friend Eigen::SparseMatrix<float> derivative0(Mesh2D& mesh,bool dual);
    friend Eigen::SparseMatrix<float> derivative1(Mesh2D& mesh,bool dual);

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& getBasisField();
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& getBasisCoefficients();

    void integrate();
private:
    unsigned int width;
    unsigned int height;
    unsigned int resX;
    unsigned int resY;

    void voxelize(unsigned char* pixels);
    bool checkVoxel(unsigned char* offs);

    bool* voxelMap;

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenValues;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenVectors;

public:
    Eigen::VectorXf velocityField;
    Eigen::SparseMatrix<float> b2;
    Eigen::SparseMatrix<float> b1;

    std::set<std::tuple<unsigned int,unsigned int,unsigned int>> faces;
    std::set<std::tuple<unsigned int,unsigned int>> edges;
    std::set<unsigned int> points;
    std::map<unsigned int,Vertex2D> vertex;
};

#endif // MESH2D_H
