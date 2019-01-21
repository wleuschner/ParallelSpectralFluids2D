#include "dec.h"
#include <iostream>


Eigen::SparseMatrix<float> hodge2(DECMesh2D& mesh,float area,bool dual)
{
    Eigen::SparseMatrix<float> h,d;
    if(dual)
    {
        d = derivative0(mesh,false);
    }
    float primalV = area;
    if(dual)
    {
        h.resize(mesh.getNumPoints(),mesh.getNumPoints());
        for(unsigned int i=0;i<mesh.getNumPoints();i++)
        {
            unsigned int nVerts=0;
            for(Eigen::SparseMatrix<float>::InnerIterator it(d,i);it;++it)
            {
                nVerts++;
            }
            if(nVerts==2)
            {
                h.insert(i,i)=1.0f/(primalV/4);
            }
            else if(nVerts==3)
            {
                h.insert(i,i)=1.0f/(primalV/2);
            }
            else if(nVerts==4)
            {
                h.insert(i,i)=1.0f/(primalV);
            }
        }
    }
    else
    {
        h.resize(mesh.getNumFaces(),mesh.getNumFaces());
        for(unsigned int i=0;i<mesh.getNumFaces();i++)
        {
            h.insert(i,i)=1.0f/primalV;
        }
    }

    return h;
}

Eigen::SparseMatrix<float> hodge1(DECMesh2D& mesh,float area,bool dual)
{
    Eigen::SparseMatrix<float> h,d;
    if(dual)
    {
        d = derivative1(mesh,true);
    }
    else
    {
        d = derivative1(mesh);
    }
    h.resize(mesh.getNumEdges(),mesh.getNumEdges());
    for(EdgeIterator it=mesh.getEdgeIteratorBegin();it!=mesh.getEdgeIteratorEnd();it++)
    {
        unsigned int i = mesh.getEdgeIndex(*it);
        float primal=area;
        float dual=0;
        for(Eigen::SparseMatrix<float>::InnerIterator it(d,i);it;++it)
        {
            dual+=primal/2;
        }
        if(dual)
        {
            h.insert(i,i)=primal/dual;
        }
        else
        {
            h.insert(i,i)=dual/primal;
        }
    }
    return h;
}

Eigen::SparseMatrix<float> hodge0(DECMesh2D& mesh,float area,bool dual)
{
    Eigen::SparseMatrix<float> h;
    float primalV = 1.0f;
    float dualV = area;
    if(dual)
    {
        h.resize(mesh.getNumFaces(),mesh.getNumFaces());
        for(unsigned int i=0;i<mesh.getNumFaces();i++)
        {
            h.insert(i,i)=primalV/dualV;
        }
    }
    else
    {
        h.resize(mesh.getNumPoints(),mesh.getNumPoints());
        for(unsigned int i=0;i<mesh.getNumPoints();i++)
        {
            h.insert(i,i)=dualV/primalV;
        }
    }

    return h;
}

Eigen::SparseMatrix<float> derivative1(DECMesh2D& mesh,bool dual)
{
    if(dual)
    {
        return derivative0(mesh).transpose();
    }
    Eigen::SparseMatrix<float> d;
    d.resize(mesh.getNumEdges(),mesh.getNumFaces());
    for(FaceIterator it = mesh.getFaceIteratorBegin();it!=mesh.getFaceIteratorEnd();it++)
    {
        unsigned int v1 = std::get<0>(*it);
        unsigned int v2 = std::get<1>(*it);
        unsigned int v3 = std::get<2>(*it);
        unsigned int v4 = std::get<3>(*it);

        d.insert(mesh.getEdgeIndex(v1,v2),mesh.getFaceIndex(*it)) = mesh.getEdgeSignum(v1,v2);
        d.insert(mesh.getEdgeIndex(v2,v3),mesh.getFaceIndex(*it)) = mesh.getEdgeSignum(v2,v3);
        d.insert(mesh.getEdgeIndex(v3,v4),mesh.getFaceIndex(*it)) = mesh.getEdgeSignum(v3,v4);
        d.insert(mesh.getEdgeIndex(v4,v1),mesh.getFaceIndex(*it)) = mesh.getEdgeSignum(v4,v1);
    }
    return d.transpose();
}

Eigen::SparseMatrix<float> derivative0(DECMesh2D& mesh,bool dual)
{
    if(dual)
    {
        return derivative1(mesh,false).transpose();
    }
    Eigen::SparseMatrix<float> d;
    d.resize(mesh.getNumPoints(),mesh.getNumEdges());
    for(EdgeIterator it = mesh.getEdgeIteratorBegin();it!=mesh.getEdgeIteratorEnd();it++)
    {
        unsigned int v1 = std::get<0>(*it);
        unsigned int v2 = std::get<1>(*it);
        d.insert(mesh.getPointIndex(v1),mesh.getEdgeIndex(*it)) = 1.0f;
        d.insert(mesh.getPointIndex(v2),mesh.getEdgeIndex(*it)) = -1.0f;
    }
    return d.transpose();
}
