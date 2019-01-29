#ifndef DECMESH2D_H
#define DECMESH2D_H
#include<set>
#include<vector>
#include"face2d.h"
#include"edge2d.h"
#include"vertex2d.h"

typedef std::vector<Vertex2D>::iterator PointIterator;
typedef std::vector<Edge2D>::iterator EdgeIterator;
typedef std::vector<Face2D>::iterator FaceIterator;

class DECMesh2D
{
public:
    DECMesh2D();
    DECMesh2D(unsigned int resolution,unsigned int voxelSize);

    void addPoint(const Vertex2D& v);
    void addEdge(const Edge2D& e);
    void addFace(const Face2D& f);

    PointIterator& getPointIteratorBegin();
    EdgeIterator& getEdgeIteratorBegin();
    FaceIterator& getFaceIteratorBegin();

    PointIterator& getPointIteratorEnd();
    EdgeIterator& getEdgeIteratorEnd();
    FaceIterator& getFaceIteratorEnd();

    unsigned int getPointIndex(const Vertex2D& v);
    unsigned int getEdgeIndex(const Edge2D& e);
    unsigned int getFaceIndex(const Face2D& f);

    bool isPointInside(const glm::vec2& point);

    Vertex2D getPoint(unsigned int id);
    Edge2D getEdge(unsigned int id);
    Face2D getFace(unsigned int id);

    int getPointSignum(const Vertex2D& v);
    int getEdgeSignum(const Edge2D& e);
    int getFaceSignum(const Face2D& f);

    int getEdgeSignum(unsigned int id,unsigned int v1,unsigned int v2);

    std::tuple<unsigned int,unsigned int,
               unsigned int,unsigned int> getIntepolationIndices(glm::vec2 coords);

    unsigned int getNumPoints();
    unsigned int getNumEdges();
    unsigned int getNumFaces();
private:
    unsigned int voxelSize;
    unsigned int resolution;

    PointIterator pointsBegin;
    EdgeIterator edgesBegin;
    FaceIterator facesBegin;

    PointIterator pointsEnd;
    EdgeIterator edgesEnd;
    FaceIterator facesEnd;

    std::vector<Vertex2D> points;
    std::vector<Edge2D> edges;
    std::vector<Face2D> faces;
};

#endif // DECMESH2D_H
