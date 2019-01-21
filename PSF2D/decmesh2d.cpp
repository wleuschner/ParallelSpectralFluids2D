#include "decmesh2d.h"
#include <tuple>

DECMesh2D::DECMesh2D()
{

}

void DECMesh2D::addPoint(unsigned int v1)
{
    if((points.find(v1))==points.end())
    {
        points.emplace(v1);
    }
}

void DECMesh2D::addEdge(unsigned int v1,unsigned int v2)
{
    if((edges.find(std::make_tuple(v1,v2)))==edges.end()&&
       (edges.find(std::make_tuple(v2,v1)))==edges.end())
    {
        edges.emplace(std::make_tuple(v1,v2));
        addPoint(v1);
        addPoint(v2);
    }
}

void DECMesh2D::addFace(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4)
{
    faces.emplace(std::make_tuple(v1,v2,v3,v4));
    addEdge(v1,v2);
    addEdge(v2,v3);
    addEdge(v3,v4);
    addEdge(v4,v1);
}

PointIterator& DECMesh2D::getPointIteratorBegin()
{
    pointsBegin = points.begin();
    return pointsBegin;
}

EdgeIterator& DECMesh2D::getEdgeIteratorBegin()
{
    edgesBegin = edges.begin();
    return edgesBegin;
}

FaceIterator& DECMesh2D::getFaceIteratorBegin()
{
    facesBegin = faces.begin();
    return facesBegin;
}

PointIterator& DECMesh2D::getPointIteratorEnd()
{
    pointsEnd = points.end();
    return pointsEnd;
}

EdgeIterator& DECMesh2D::getEdgeIteratorEnd()
{
    edgesEnd = edges.end();
    return edgesEnd;
}

FaceIterator& DECMesh2D::getFaceIteratorEnd()
{
    facesEnd = faces.end();
    return facesEnd;
}

unsigned int DECMesh2D::getPointIndex(unsigned int v1)
{
    std::set<unsigned int>::iterator it;
    it=points.find(v1);
    return std::distance(points.begin(),it);
}

unsigned int DECMesh2D::getEdgeIndex(unsigned int v1,unsigned int v2)
{
    std::set<std::tuple<unsigned int,unsigned int>>::iterator it;
    if((it=edges.find(std::make_tuple(v1,v2)))==edges.end())
    {
        it=edges.find(std::make_tuple(v2,v1));
    }
    return std::distance(edges.begin(),it);
}

unsigned int DECMesh2D::getFaceIndex(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4)
{
    std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>>::iterator it;
    it=faces.find(std::make_tuple(v1,v2,v3,v4));
    return std::distance(faces.begin(),it);
}

unsigned int DECMesh2D::getEdgeIndex(std::tuple<unsigned int,unsigned int> e)
{
    std::set<std::tuple<unsigned int,unsigned int>>::iterator it;
    if((it=edges.find(e))==edges.end())
    {
        it=edges.find(std::make_tuple(std::get<1>(e),std::get<0>(e)));
    }
    return std::distance(edges.begin(),it);
}

unsigned int DECMesh2D::getFaceIndex(std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> f)
{
    std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>>::iterator it;
    it=faces.find(f);
    return std::distance(faces.begin(),it);
}

unsigned int DECMesh2D::getPoint(unsigned int v1)
{
    std::set<unsigned int>::iterator it;
    it=points.find(v1);
    return *it;
}

std::tuple<unsigned int,unsigned int> DECMesh2D::getEdge(unsigned int v1,unsigned int v2)
{
    std::set<std::tuple<unsigned int,unsigned int>>::iterator it;
    if((it=edges.find(std::make_tuple(v1,v2)))==edges.end())
    {
        it=edges.find(std::make_tuple(v2,v1));
    }
    return *it;
}

std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> DECMesh2D::getFace(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4)
{
    std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>>::iterator it;
    it=faces.find(std::make_tuple(v1,v2,v3,v4));
    return *it;
}

int DECMesh2D::getPointSignum(unsigned int v1)
{
    return 1;
}

int DECMesh2D::getEdgeSignum(unsigned int v1,unsigned int v2)
{
    std::set<std::tuple<unsigned int,unsigned int>>::iterator it;
    if(edges.find(std::make_tuple(v1,v2))==edges.end())
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

int DECMesh2D::getFaceSignum(unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4)
{
    return 1;
}

int DECMesh2D::getEdgeSignum(std::tuple<unsigned int,unsigned int> e)
{
    std::set<std::tuple<unsigned int,unsigned int>>::iterator it;
    if(edges.find(e)==edges.end())
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

int getFaceSignum(std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> f)
{
    return 1;
}

unsigned int DECMesh2D::getNumPoints()
{
    return static_cast<unsigned int>(points.size());
}

unsigned int DECMesh2D::getNumEdges()
{
    return static_cast<unsigned int>(edges.size());
}

unsigned int DECMesh2D::getNumFaces()
{
    return static_cast<unsigned int>(faces.size());
}
