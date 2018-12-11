#include "dec.h"
#include <iostream>


Eigen::SparseMatrix<float> hodge2(Mesh2D& mesh,bool dual)
{
    Eigen::SparseMatrix<float> h;
    if(dual)
    {
        h.resize(mesh.points.size(),mesh.points.size());
        for(unsigned int i=0;i<mesh.points.size();i++)
        {
            h.insert(i,i)=64.0f;
        }
    }
    else
    {
        h.resize(mesh.faces.size(),mesh.faces.size());
        for(unsigned int i=0;i<mesh.points.size();i++)
        {
            h.insert(i,i)=1.0/64.0f;
        }
    }

    return h;
}

Eigen::SparseMatrix<float> hodge1(Mesh2D& mesh,bool inverse)
{
    Eigen::SparseMatrix<float> h,d;
    if(inverse)
    {
        d = derivative1(mesh,true);
    }
    else
    {
        d = derivative1(mesh);
    }
    h.resize(mesh.edges.size(),mesh.edges.size());
    for(unsigned int i=0;i<mesh.edges.size();i++)
    {
        float dual=0;
        for(Eigen::SparseMatrix<float>::InnerIterator it(d,i);it;++it)
        {
            dual+=16.0f;
        }
        std::cout<<"DUAL: "<<dual<<std::endl;
        if(inverse)
        {
            h.insert(i,i)=32.0f/dual;
        }
        else
        {
            h.insert(i,i)=dual/32.0f;
        }
    }
    return h;
}

Eigen::SparseMatrix<float> hodge0(Mesh2D& mesh,bool inverse)
{
    Eigen::SparseMatrix<float> h;
    if(inverse)
    {
        h.resize(mesh.faces.size(),mesh.faces.size());
    }
    else
    {
        h.resize(mesh.points.size(),mesh.points.size());
    }
    return h;
}

//TODO implement derivative
Eigen::SparseMatrix<float> derivative1(Mesh2D& mesh,bool dual)
{
    if(dual)
    {
        return derivative0(mesh).transpose();
    }
    Eigen::SparseMatrix<float> d;
    d.resize(mesh.edges.size(),mesh.faces.size());
    std::cout<<mesh.edges.size()<<" "<<mesh.faces.size();
    for(auto it = mesh.faces.begin();it!=mesh.faces.end();it++)
    {
        std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> f=*it;
        unsigned int v1 = std::get<0>(f);
        unsigned int v2 = std::get<1>(f);
        unsigned int v3 = std::get<2>(f);
        unsigned int v4 = std::get<3>(f);
        std::set<std::tuple<unsigned int,unsigned int>>::iterator eit;
        if((eit=mesh.edges.find(std::make_tuple(v1,v2)))!=mesh.edges.end())
        {
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=1.0f;
        }
        else
        {
            eit=mesh.edges.find(std::make_tuple(v2,v1));
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=-1.0f;
        }

        if((eit=mesh.edges.find(std::make_tuple(v2,v3)))!=mesh.edges.end())
        {
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=1.0f;
        }
        else
        {
            eit=mesh.edges.find(std::make_tuple(v3,v2));
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=-1.0f;
        }

        if((eit=mesh.edges.find(std::make_tuple(v3,v4)))!=mesh.edges.end())
        {
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=1.0f;
        }
        else
        {
            eit=mesh.edges.find(std::make_tuple(v4,v3));
            std::cout<<std::distance(mesh.edges.begin(),eit)<<" "<<std::distance(mesh.faces.begin(),it)<<std::endl;
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=-1.0f;
        }
        if((eit=mesh.edges.find(std::make_tuple(v4,v1)))!=mesh.edges.end())
        {
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=1.0f;
        }
        else
        {
            eit=mesh.edges.find(std::make_tuple(v1,v4));
            std::cout<<std::distance(mesh.edges.begin(),eit)<<" "<<std::distance(mesh.faces.begin(),it)<<std::endl;
            d.insert(std::distance(mesh.edges.begin(),eit),std::distance(mesh.faces.begin(),it))=-1.0f;
        }
    }
    return d.transpose();
}

Eigen::SparseMatrix<float> derivative0(Mesh2D& mesh,bool dual)
{
    if(dual)
    {
        return derivative1(mesh,false).transpose();
    }
    Eigen::SparseMatrix<float> d;
    d.resize(mesh.points.size(),mesh.edges.size());
    for(auto it = mesh.edges.begin();it!=mesh.edges.end();it++)
    {
        std::tuple<unsigned int,unsigned int> f=*it;
        unsigned int v1 = std::get<0>(f);
        unsigned int v2 = std::get<1>(f);
        std::set<unsigned int>::iterator pit1,pit2;
        if((pit1=mesh.points.find(v1))!=mesh.points.end())
        {
            d.insert(std::distance(mesh.points.begin(),pit1),std::distance(mesh.edges.begin(),it))=1;
        }
        if((pit2=mesh.points.find(v2))!=mesh.points.end())
        {
            d.insert(std::distance(mesh.points.begin(),pit2),std::distance(mesh.edges.begin(),it))=-1;
        }
    }
    return d.transpose();
}
