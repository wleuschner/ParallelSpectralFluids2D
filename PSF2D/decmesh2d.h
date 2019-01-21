#ifndef DECMESH2D_H
#define DECMESH2D_H
#include<set>

typedef std::set<unsigned int>::iterator PointIterator;
typedef std::set<std::tuple<unsigned int,unsigned int>>::iterator EdgeIterator;
typedef std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>>::iterator FaceIterator;

class DECMesh2D
{
public:
    DECMesh2D();

    void addPoint(unsigned int v1);
    void addEdge(unsigned int v1,unsigned int v2);
    void addFace(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4);

    PointIterator& getPointIteratorBegin();
    EdgeIterator& getEdgeIteratorBegin();
    FaceIterator& getFaceIteratorBegin();

    PointIterator& getPointIteratorEnd();
    EdgeIterator& getEdgeIteratorEnd();
    FaceIterator& getFaceIteratorEnd();

    unsigned int getPointIndex(unsigned int v1);
    unsigned int getEdgeIndex(unsigned int v1,unsigned int v2);
    unsigned int getFaceIndex(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4);

    unsigned int getEdgeIndex(std::tuple<unsigned int,unsigned int> e);
    unsigned int getFaceIndex(std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> f);

    unsigned int getPoint(unsigned int v1);
    std::tuple<unsigned int,unsigned int> getEdge(unsigned int v1,unsigned int v2);
    std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> getFace(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4);

    int getPointSignum(unsigned int v1);
    int getEdgeSignum(unsigned int v1,unsigned int v2);
    int getFaceSignum(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4);

    int getEdgeSignum(std::tuple<unsigned int,unsigned int> e);
    int getFaceSignum(std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> f);

    unsigned int getNumPoints();
    unsigned int getNumEdges();
    unsigned int getNumFaces();
private:
    PointIterator pointsBegin;
    EdgeIterator edgesBegin;
    FaceIterator facesBegin;

    PointIterator pointsEnd;
    EdgeIterator edgesEnd;
    FaceIterator facesEnd;

    std::set<unsigned int> points;
    std::set<std::tuple<unsigned int,unsigned int>> edges;
    std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>> faces;
};

#endif // DECMESH2D_H
