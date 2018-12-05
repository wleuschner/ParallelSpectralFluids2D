#include "dec.h"
#include <iostream>


Eigen::SparseMatrix<float> hodge2(Mesh2D& mesh,bool dual)
{
    Eigen::SparseMatrix<float> ret;
    unsigned int nFaces;
    if(dual)
    {
        ret.resize(mesh.dual->faces.size(),mesh.dual->faces.size());
        nFaces = mesh.dual->faces.size();
    }
    else
    {
        ret.resize(mesh.faces.size(),mesh.faces.size());
        nFaces = mesh.faces.size();
    }
    for(unsigned int k=0;k<nFaces;++k)
    {
        if(dual)
        {
            if(mesh.dual->faces[k].area==0.0f)
            {
                continue;
            }
            ret.insert(k,k) = mesh.vertex[k].area/mesh.dual->faces[k].area;
        }
        else
        {
            if(mesh.faces[k].area==0.0f)
            {
                continue;
            }
            ret.insert(k,k) = mesh.dual->vertex[k].area/mesh.faces[k].area;
        }
    }
    std::cout<<"HODGE2 DIM: "<<ret.cols()<<" "<<ret.rows()<<std::endl;
    return ret;
}

Eigen::SparseMatrix<float> hodge1(Mesh2D& mesh,bool inverse)
{
    Eigen::SparseMatrix<float> ret;
    ret.resize(mesh.edges.size(),mesh.edges.size());
    for(unsigned int k=0;k<mesh.edges.size();++k)
    {
        if(mesh.edges[k].area==0.0f || mesh.dual->edges[k].area==0.0f)
        {
            continue;
        }
        else if(inverse)
        {
            ret.insert(k,k) = mesh.edges[k].area/mesh.dual->edges[k].area;
            std::cout<<"("<<k<<","<<k<<")"<<"="<<mesh.edges[k].area/mesh.dual->edges[k].area<<std::endl;
        }
        else
        {
            ret.insert(k,k) = mesh.dual->edges[k].area/mesh.edges[k].area;
            std::cout<<"("<<k<<","<<k<<")"<<"="<<mesh.dual->edges[k].area/mesh.edges[k].area<<std::endl;
        }
    }
    return ret;
}

Eigen::SparseMatrix<float> hodge0(Mesh2D& mesh,bool inverse)
{
    Eigen::SparseMatrix<float> ret;
    ret.resize(mesh.vertex.size(),mesh.vertex.size());
    for(unsigned int k=0;k<mesh.vertex.size();++k)
    {
        ret.insert(k,k) = mesh.dual->faces[k].area/mesh.vertex[k].area;
    }
    return ret;
}

//TODO implement derivative
Eigen::SparseMatrix<float> derivative0(Mesh2D& mesh,bool dual)
{
    Eigen::SparseMatrix<float> d;
    unsigned int nEdges;
    if(dual)
    {
        d.resize(mesh.dual->edges.size(),mesh.dual->vertex.size());
        nEdges = mesh.dual->edges.size();
    }
    else
    {
        d.resize(mesh.edges.size(),mesh.vertex.size());
        nEdges = mesh.edges.size();
    }

    for(unsigned int e=0;e<nEdges;e++)
    {
        Edge2D edge;
        if(dual)
        {
            edge = mesh.dual->edges[e];
        }
        else
        {
            edge = mesh.edges[e];
        }
        if(edge.area!=0.0f)
        {
            if(edge.v1!=-1)
            {
                d.insert(e,edge.v1) = 1;
            }
            if(edge.v2!=-1)
            {
                d.insert(e,edge.v2) = 1;
            }
        }
    }
    return d;
}

Eigen::SparseMatrix<float> derivative1(Mesh2D& mesh,bool dual)
{
    Eigen::SparseMatrix<float> d;
    unsigned int nFaces;
    if(dual)
    {
        d.resize(mesh.dual->faces.size(),mesh.dual->edges.size());
        nFaces = mesh.dual->faces.size();
    }
    else
    {
        d.resize(mesh.faces.size(),mesh.edges.size());
        nFaces = mesh.faces.size();
    }
    for(unsigned int f=0;f<nFaces;f++)
    {
        Face2D face;
        if(dual)
        {
            face = mesh.dual->faces[f];
        }
        else
        {
            face = mesh.faces[f];
        }
        if(face.area!=0.0f)
        {
            d.insert(f,std::abs(face.e1)) = face.e1>0?1:-1;
            d.insert(f,std::abs(face.e2)) = face.e2>0?1:-1;
            d.insert(f,std::abs(face.e3)) = face.e3>0?1:-1;
            d.insert(f,std::abs(face.e4)) = face.e4>0?1:-1;
        }
    }
    return d;
}
