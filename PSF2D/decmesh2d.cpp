#include "decmesh2d.h"
#include <tuple>

DECMesh2D::DECMesh2D()
{
}

DECMesh2D::DECMesh2D(unsigned int resolution)
{
    this->resolution = resolution;
    faces.resize(resolution*resolution);
    edges.resize(2*(resolution+1)*(resolution+1)-2*(resolution+1));
    points.resize((resolution+1)*(resolution+1));
}

void DECMesh2D::addPoint(const Vertex2D& v)
{
    if(points[v.id].inside==GridState::UNINITIALIZED)
    {
        points[v.id] = v;
    }
    else if(points[v.id].inside==GridState::OUTSIDE&&
            v.inside==GridState::INSIDE)
    {
        points[v.id].inside = v.inside;
    }
}

void DECMesh2D::addEdge(const Edge2D& e)
{
    if(edges[e.id].inside==GridState::UNINITIALIZED)
    {
        edges[e.id] = e;
        addPoint(Vertex2D(e.v1,e.inside));
        addPoint(Vertex2D(e.v2,e.inside));
    }
    else if(edges[e.id].inside==GridState::OUTSIDE&&
            e.inside==GridState::INSIDE)
    {
        edges[e.id].inside = e.inside;
    }
}

void DECMesh2D::addFace(const Face2D& f)
{
    if(faces[f.id].inside==GridState::UNINITIALIZED)
    {
        faces[f.id] = f;
        unsigned int yOfs = f.id/resolution;
        unsigned int xOfs = f.id%resolution;

        unsigned int eid1 = yOfs*(resolution+resolution+1)+resolution+xOfs;
        unsigned int eid2 = (yOfs+1)*(resolution+resolution+1)+xOfs;
        unsigned int eid3 = yOfs*(resolution+resolution+1)+resolution+xOfs+1;
        unsigned int eid4 = yOfs*(resolution+resolution+1)+xOfs;

        faces[f.id].e1 = eid1;
        faces[f.id].e2 = eid2;
        faces[f.id].e3 = eid3;
        faces[f.id].e4 = eid4;

        addEdge(Edge2D(eid1,f.v1,f.v2,f.inside));
        addEdge(Edge2D(eid2,f.v2,f.v3,f.inside));
        addEdge(Edge2D(eid3,f.v3,f.v4,f.inside));
        addEdge(Edge2D(eid4,f.v4,f.v1,f.inside));
    }
    else if(faces[f.id].inside==GridState::OUTSIDE&&
            f.inside==GridState::INSIDE)
    {
        faces[f.id].inside = f.inside;
    }
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

unsigned int DECMesh2D::getPointIndex(const Vertex2D& v)
{
    return v.id;
}

unsigned int DECMesh2D::getEdgeIndex(const Edge2D& e)
{
    return e.id;
}

unsigned int DECMesh2D::getFaceIndex(const Face2D& f)
{
    return f.id;
}

Vertex2D DECMesh2D::getPoint(unsigned int id)
{
    return points[id];
}

Edge2D DECMesh2D::getEdge(unsigned int id)
{
    return edges[id];
}

Face2D DECMesh2D::getFace(unsigned int id)
{
    return faces[id];
}

int DECMesh2D::getPointSignum(const Vertex2D& v)
{
    return 1;
}

int DECMesh2D::getEdgeSignum(const Edge2D& e)
{
    if(edges[e.id].v1==e.v1 && edges[e.id].v2==e.v2)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

int DECMesh2D::getFaceSignum(const Face2D& f)
{
    return 1;
}

int DECMesh2D::getEdgeSignum(unsigned int id,unsigned int v1,unsigned int v2)
{
    if(edges[id].v1==v1 && edges[id].v2==v2)
    {
        return 1;
    }
    else
    {
        return -1;
    }
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
