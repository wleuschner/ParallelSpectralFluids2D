#include "dec.h"
#include <iostream>


Eigen::SparseMatrix<double> hodge2(DECMesh2D& mesh,double area,bool dual)
{
    Eigen::SparseMatrix<double> h,b0,b1;

    if(dual)
    {
        b0 = derivative0(mesh);
        b1 = derivative1(mesh);
        h.resize(mesh.getNumPoints(),mesh.getNumPoints());
        for(PointIterator pit=mesh.getPointIteratorBegin();pit!=mesh.getPointIteratorEnd();pit++)
        {
            unsigned int i=mesh.getPointIndex(*pit);
            unsigned int nFaces=0;
            for(Eigen::SparseMatrix<double>::InnerIterator it(b0,i);it;++it)
            {
                for(Eigen::SparseMatrix<double>::InnerIterator fit(b1,it.row());fit;++fit)
                {
                    nFaces++;
                }
            }
            if(nFaces==2)
            {
                h.insert(i,i)=4.0;
            }
            else if(nFaces==4)
            {
                h.insert(i,i)=2.0;
            }
            else if(nFaces==6)
            {
                h.insert(i,i)=1.3333333333333;
            }
            else if(nFaces==8)
            {
                h.insert(i,i)=1.0;
            }
        }
    }
    else
    {
        h.resize(mesh.getNumFaces(),mesh.getNumFaces());
        for(unsigned int i=0;i<mesh.getNumFaces();i++)
        {
            h.insert(i,i)=1.0/1.0;
        }
    }

    return h;
}

Eigen::SparseMatrix<double> hodge1(DECMesh2D& mesh,double area,bool dual)
{
    Eigen::SparseMatrix<double> h,d;

    d = derivative1(mesh);

    h.resize(mesh.getNumEdges(),mesh.getNumEdges());
    h.setIdentity();
    return h;

    for(EdgeIterator eit=mesh.getEdgeIteratorBegin();eit!=mesh.getEdgeIteratorEnd();eit++)
    {
        unsigned int i = mesh.getEdgeIndex(*eit);
        unsigned int dualEdges=0;
        for(Eigen::SparseMatrix<double>::InnerIterator it(d,i);it;++it)
        {
            dualEdges++;
        }
        if(dual)
        {
            if(dualEdges==1)
            {
                h.insert(i,i)=2.0;
            }
            else if(dualEdges==2)
            {
                h.insert(i,i)=1.0;
            }
        }
        else
        {
            if(dualEdges==1)
            {
                h.insert(i,i)=0.5;
            }
            else if(dualEdges==2)
            {
                h.insert(i,i)=1.0;
            }
        }
    }
    return h;
}

Eigen::SparseMatrix<double> hodge0(DECMesh2D& mesh,double area,bool dual)
{
    Eigen::SparseMatrix<double> h;
    double primalV = 1.0f;
    double dualV = area;
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

Eigen::SparseMatrix<double> derivative1(DECMesh2D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> d;
    if(dual)
    {
        return derivative0(mesh,false).transpose();
    }
    d.resize(mesh.getNumEdges(),mesh.getNumFaces());
    for(FaceIterator it = mesh.getFaceIteratorBegin();it!=mesh.getFaceIteratorEnd();it++)
    {
        unsigned int v1 = std::get<0>(*it);
        unsigned int v2 = std::get<1>(*it);
        unsigned int v3 = std::get<2>(*it);
        unsigned int v4 = std::get<3>(*it);
        unsigned int fidx = mesh.getFaceIndex(*it);

        d.insert(mesh.getEdgeIndex(v1,v2),fidx) = mesh.getEdgeSignum(v1,v2);
        d.insert(mesh.getEdgeIndex(v2,v3),fidx) = mesh.getEdgeSignum(v2,v3);
        d.insert(mesh.getEdgeIndex(v3,v4),fidx) = mesh.getEdgeSignum(v3,v4);
        d.insert(mesh.getEdgeIndex(v4,v1),fidx) = mesh.getEdgeSignum(v4,v1);
    }
    return d.transpose();
}

Eigen::SparseMatrix<double> derivative0(DECMesh2D& mesh,bool dual)
{
    if(dual)
    {
        return derivative1(mesh,false).transpose();
    }
    Eigen::SparseMatrix<double> d;
    d.resize(mesh.getNumPoints(),mesh.getNumEdges());
    for(EdgeIterator it = mesh.getEdgeIteratorBegin();it!=mesh.getEdgeIteratorEnd();it++)
    {
        unsigned int v1 = std::get<0>(*it);
        unsigned int v2 = std::get<1>(*it);
        unsigned int eidx = mesh.getEdgeIndex(*it);
        d.insert(mesh.getPointIndex(v1),eidx) = -1.0f;
        d.insert(mesh.getPointIndex(v2),eidx) = 1.0f;
    }
    return d.transpose();
}
