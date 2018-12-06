#include "mesh2d.h"
#include <iostream>
#include <eigen3/unsupported/Eigen/ArpackSupport>
#include "dec.h"
typedef Eigen::SparseMatrix<float> SparseMat;
typedef Eigen::SimplicialLDLT<SparseMat> SparseChol;
typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver <SparseMat, SparseChol> Arpack;

Mesh2D::Mesh2D()
{
}

Mesh2D::Mesh2D(unsigned int w,unsigned int h,unsigned int resX,unsigned int resY,unsigned char* data) : Mesh2D()
{
    this->width = w;
    this->height = h;
    this->resX = resX;
    this->resY = resY;
    voxelize(data);
    buildLaplace();
    for(unsigned int y=0;y<height/resY;y++)
    {
        for(unsigned int x=0;x<width/resX;x++)
        {
            std::cout<<voxelMap[y*width/resX+x];
        }
        std::cout<<std::endl;
    }
}

Mesh2D::~Mesh2D()
{
}

unsigned int Mesh2D::getNumVertices()
{
    return vertex.size();
}

unsigned int Mesh2D::getNumEdges()
{
    return edges.size();
}

unsigned int Mesh2D::getNumFaces()
{
    return faces.size();
}

bool* Mesh2D::getVoxelMap()
{
    return voxelMap;
}

unsigned int Mesh2D::getWidth()
{
    return width;
}

unsigned int Mesh2D::getHeight()
{
    return height;
}

unsigned int Mesh2D::getResX()
{
    return resX;
}

unsigned int Mesh2D::getResY()
{
    return resY;
}

Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& Mesh2D::getBasisField()
{
    return eigenVectors;
}

Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& Mesh2D::getBasisCoefficients()
{
    return eigenValues;
}

Eigen::VectorXf Mesh2D::getVelocityField()
{
    return velocityField;
}

void Mesh2D::voxelize(unsigned char* pixels)
{
    voxelMap = new bool[width/resX*height/resY];
/*    vertex.resize((height/resY+1)*(width/resX+1));
    edges.resize(2*(height/resY+2)*(width/resX+1)-(height/resY+2)-(width/resX+1));
    faces.resize(height/resY*width/resX);*/

    unsigned int nFaces = height/resY*width/resX;
    unsigned int nVerts = (height/resY+1)*(width/resX+1);
    unsigned int nEdges = 2*(height/resY+2)*(width/resX+1)-(height/resY+2)-(width/resX+1);
    b2.reserve(nEdges*nFaces);
    b1.reserve(nVerts*nEdges);
    //faces.resize(width/w*height/h,width/w*height/h);
    unsigned int iv1,iv2,iv3,iv4;
    int e1,e2,e3,e4;
    int f;
    int y;


    for(int y=0;y<height/resY;y++)
    {
        int x;
        for(int x=0;x<width/resX;x++)
        {
            unsigned char* start = &pixels[y*width*resY+x*resX];
            iv1 = y*(width/resX)+x; //UP LEFT Vert
            iv2 = y*(width/resX)+x+1; //UP RIGHT Vert
            iv3 = (y+1)*(width/resX)+x; //DOWN LEFT vert
            iv4 = (y+1)*(width/resX)+x+1; //DOWN RIGHT vert*/

            f = y*(width/resX)+x; //FACE
            if(checkVoxel(start))
            {
                voxelMap[f] = true;
                Vertex2D v1,v2,v3,v4;
                v1.pos.x = x*resX;
                v1.pos.y = y*resY;
                v2.pos.x = (x+1)*resX;
                v2.pos.y = y*resY;
                v3.pos.x = x*resX;
                v3.pos.y = (y+1)*resY;
                v4.pos.x = (x+1)*resX;
                v4.pos.y = (y+1)*resY;
                vertex[iv1] = v1;
                vertex[iv2] = v2;
                vertex[iv3] = v3;
                vertex[iv4] = v4;
                std::set<std::tuple<unsigned int,unsigned int,unsigned int>>::const_iterator itf1=faces.emplace(std::make_tuple(iv1,iv3,iv2)).first;
                std::set<std::tuple<unsigned int,unsigned int,unsigned int>>::const_iterator itf2=faces.emplace(std::make_tuple(iv3,iv4,iv2)).first;
                std::set<std::tuple<unsigned int,unsigned int>>::const_iterator it;
                std::set<unsigned int>::const_iterator pit1,pit2;

                //First Triangle
                if((it=edges.find(std::make_tuple(iv1,iv3)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    pit1=points.find(iv1);
                    pit2=points.find(iv3);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else if((it=edges.find(std::make_tuple(iv3,iv1)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                    pit1=points.find(iv3);
                    pit2=points.find(iv1);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;

                }
                else
                {
                    it = edges.emplace(std::make_tuple(iv1,iv3)).first;
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    if((pit1=points.find(iv1))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    }
                    else
                    {
                        pit1 = points.emplace(iv1).first;
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        vertex[iv1] = v1;
                    }
                    if((pit2=points.find(iv3))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        pit2 = points.emplace(iv3).first;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        vertex[iv3] = v3;
                    }
                }

                if((it=edges.find(std::make_tuple(iv2,iv1)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    pit1=points.find(iv2);
                    pit2=points.find(iv1);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else if((it=edges.find(std::make_tuple(iv1,iv2)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                    pit1=points.find(iv1);
                    pit2=points.find(iv2);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else
                {
                    it = edges.emplace(std::make_tuple(iv2,iv1)).first;
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    if((pit1=points.find(iv2))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    }
                    else
                    {
                        pit1 = points.emplace(iv2).first;
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        vertex[iv2] = v2;
                    }
                    if((pit2=points.find(iv1))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        pit2 = points.emplace(iv1).first;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        vertex[iv1] = v1;
                    }
                }

                if((it=edges.find(std::make_tuple(iv3,iv2)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    pit1=points.find(iv3);
                    pit2=points.find(iv2);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else if((it=edges.find(std::make_tuple(iv2,iv3)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                    pit1=points.find(iv2);
                    pit2=points.find(iv3);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else
                {
                    it = edges.emplace(std::make_tuple(iv3,iv2)).first;
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    if((pit1=points.find(iv3))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    }
                    else
                    {
                        pit1 = points.emplace(iv3).first;
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        vertex[iv3] = v3;
                    }
                    if((pit2=points.find(iv2))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        pit2 = points.emplace(iv2).first;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        vertex[iv2] = v2;
                    }
                }

                //Second Triangle
                if((it=edges.find(std::make_tuple(iv3,iv4)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    pit1=points.find(iv3);
                    pit2=points.find(iv4);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else if((it=edges.find(std::make_tuple(iv4,iv3)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                    pit1=points.find(iv4);
                    pit2=points.find(iv3);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else
                {
                    it = edges.emplace(std::make_tuple(iv3,iv4)).first;
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    if((pit1=points.find(iv3))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    }
                    else
                    {
                        pit1 = points.emplace(iv3).first;
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        vertex[iv3] = v3;
                    }
                    if((pit2=points.find(iv4))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        pit2 = points.emplace(iv4).first;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        vertex[iv4] = v4;
                    }
                }

                if((it=edges.find(std::make_tuple(iv2,iv3)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    pit1=points.find(iv2);
                    pit2=points.find(iv3);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else if((it=edges.find(std::make_tuple(iv3,iv2)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                    pit1=points.find(iv3);
                    pit2=points.find(iv2);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else
                {
                    it = edges.emplace(std::make_tuple(iv2,iv3)).first;
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    if((pit1=points.find(iv2))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    }
                    else
                    {
                        pit1 = points.emplace(iv2).first;
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        vertex[iv2] = v2;
                    }
                    if((pit2=points.find(iv3))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        pit2 = points.emplace(iv3).first;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        vertex[iv3] = v3;
                    }
                }

                if((it=edges.find(std::make_tuple(iv4,iv2)))!=edges.end())
                {
                   //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                   pit1=points.find(iv4);
                   pit2=points.find(iv2);
                   //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                   //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else if((it=edges.find(std::make_tuple(iv2,iv4)))!=edges.end())
                {
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                    pit1=points.find(iv2);
                    pit2=points.find(iv4);
                    //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                }
                else
                {
                    it = edges.emplace(std::make_tuple(iv4,iv2)).first;
                    //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                    if((pit1=points.find(iv4))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                    }
                    else
                    {
                        pit1 = points.emplace(iv4).first;
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        vertex[iv4] = v4;
                    }
                    if((pit2=points.find(iv2))!=points.end())
                    {
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        pit2 = points.emplace(iv2).first;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        vertex[iv2] = v2;
                    }
                }
            }
            else
            {
                voxelMap[f] = false;
            }
        }
    }
}

bool Mesh2D::checkVoxel(unsigned char* offs)
{
    for(unsigned int ry=0;ry<resY;ry++)
    {
        for(unsigned int rx=0;rx<resX;rx++)
        {
            if(offs[ry*width+rx]==0)
            {
                return true;
            }
        }
    }
    return false;
}

Eigen::SparseMatrix<float> Mesh2D::buildLaplace()
{
    //Build Face Matrix
    SparseMat mat,A;
    mat.resize(edges.size(),edges.size());
    mat = derivative1(*this).transpose();
    //derivative1(*this)*
    //Calculate Laplace Operator
    //mat = -1.0f*derivative1(*this)*hodge1(*this,true)*derivative1(*this).transpose()*hodge2(*this,false);
    //std::cout<<"HODGE DIMS"<<hodge1(*this,false).cols()<<" "<<hodge1(*this,false).rows()<<std::endl;
    //std::cout<<"DERIVATIVE DIMS"<<derivative1(*this,false).cols()<<" "<<derivative0(*this,false).rows()<<std::endl;

    //mat = -1.0f*derivative0(*this)*hodge2(*this,true)*derivative1(*this,true)*hodge1(*this,false);
    std::cout<<"MAT DIMS: "<<mat.cols()<<" "<<mat.rows()<<std::endl;
    //Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<float>> eigenSolver;
    //eigenSolver.compute(A);
    //mat = hodge1(*this,true);//derivative1(*this).transpose()*hodge2(*this,false);
    //std::cout<<eigenSolver.eigenvalues()<<std::endl;

    std::cout<<mat<<std::endl;

    Arpack arpack;
    arpack.compute(mat,32,"LM");
    eigenValues = arpack.eigenvalues();
    eigenVectors = arpack.eigenvectors();
    std::cout<<eigenValues.transpose()<<std::endl;
    std::cout<<eigenVectors<<std::endl;

    std::cout<<"Velocity"<<std::endl;
    std::cout.precision(5);

    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    velocityField = Eigen::VectorXf::Zero(getNumEdges());
    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        Eigen::VectorXf e = eigenVectors.col(k);
        velocityField += (eigenValues(k,0)*e).transpose();

    }
    for(unsigned int k=0;k<velocityField.rows();k++)
    {
        std::cout.width(10);
        std::cout<<velocityField(k)<<" ";
        if((k+1)%16==0)
        {
            std::cout<<std::endl;
        }
    }

    return mat;

}

void Mesh2D::integrate()
{
    std::cout<<"Advect"<<std::endl;
    /*for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        //std::cout<<(eigenValues(k,0)*eigenVectors.col(k))<<" ";
        if(k%15==0)
        {
            std::cout<<std::endl;
        }
    }*/
    float e1 = 0.0f;
    float e2 = 0.0f;
    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        e1 += eigenValues(k,0)*eigenValues(k,0);
    }
    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        e2 += eigenValues(k,0)*eigenValues(k,0);
    }
    eigenValues *= std::sqrt(e1/e2);
    std::cout<<e1<<std::endl;
}
