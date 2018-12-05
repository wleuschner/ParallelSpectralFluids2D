#include "mesh2d.h"
#include <iostream>
#include <eigen3/unsupported/Eigen/ArpackSupport>
#include "dec.h"
typedef Eigen::SparseMatrix<float> SparseMat;
typedef Eigen::SimplicialLDLT<SparseMat> SparseChol;
typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver <SparseMat, SparseChol> Arpack;

Mesh2D::Mesh2D()
{
    dual = NULL;
}

Mesh2D::Mesh2D(unsigned int w,unsigned int h,unsigned int resX,unsigned int resY,unsigned char* data) : Mesh2D()
{
    this->width = w;
    this->height = h;
    this->resX = resX;
    this->resY = resY;
    this->dual = new Mesh2D();
    voxelize(data);
    buildDual();
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
    if(dual!=NULL)
    {
        //delete dual;
        //dual = NULL;
    }
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
    vertex.resize((height/resY+1)*(width/resX+1));
    edges.resize(2*(height/resY+2)*(width/resX+1)-(height/resY+2)-(width/resX+1));
    faces.resize(height/resY*width/resX);

    //faces.resize(width/w*height/h,width/w*height/h);
    int v1,v2,v3,v4;
    int e1,e2,e3,e4;
    int f;
    int y;
    for(y=0;y<height/resY;y++)
    {
        bool first=true;
        int x;
        for(x=0;x<width/resX;x++)
        {
            unsigned char* start = &pixels[y*width*resY+x*resX];
            v1 = y*(width/resX)+x+y; //UP LEFT Vert
            v2 = y*(width/resX)+x+1+y; //UP RIGHT Vert
            v3 = (y+1)*(width/resX)+x; //DOWN LEFT vert
            v4 = (y+1)*(width/resX)+x+1; //DOWN RIGHT vert

            e1 = y*(2*width/resX)+x*2; //Edge v2->v1
            e2 = y*(2*width/resX)+x*2+1+y*2; //Edge v1->v3
            e3 = (y+1)*(2*width/resX)+x*2; //Edge v3->v4
            e4 = y*(2*width/resX)+x*2+3+y*2; //Edge v4->v2

            f = y*(width/resX)+x; //FACE
            if(checkVoxel(start))
            {
                voxelMap[f] = true;
                //faces.insert(y*width/w+x,y*width/w+x) = w*h;
            }
            else
            {
                voxelMap[f] = false;
            }
            vertex[v1].pos.x = x*resX;
            vertex[v1].pos.y = y*resY;
            vertex[v1].area = 1.0f;
            vertex[v2].pos.x = (x+1)*resX;
            vertex[v2].pos.y = y*resY;
            vertex[v2].area = 1.0f;
            vertex[v3].pos.x = x*resX;
            vertex[v3].pos.y = (y+1)*resY;
            vertex[v3].area = 1.0f;
            vertex[v4].pos.x = (x+1)*resX;
            vertex[v4].pos.y = (y+1)*resY;
            vertex[v4].area = 1.0f;

            edges[e2].faces[0] = f;
            edges[e4].faces[1] = f;
            edges[e1].faces[0] = f;
            edges[e3].faces[1] = f;

            //First Row First Column
            if(f/(width/resX)==0 && f%(width/resX)==0)
            {
                edges[e1].v1 = v2;
                edges[e1].v2 = v1;
                edges[e2].v1 = v1;
                edges[e2].v2 = v3;
                edges[e3].v1 = v3;
                edges[e3].v2 = v4;
                edges[e4].v1 = v4;
                edges[e4].v2 = v2;
                faces[f] = Face2D(e1,e2,e3,e4);
                if(checkVoxel(start))
                {
                    edges[e1].area = glm::length(vertex[edges[e1].v2].pos-vertex[edges[e1].v1].pos);
                    edges[e2].area = glm::length(vertex[edges[e2].v2].pos-vertex[edges[e2].v1].pos);
                    edges[e3].area = glm::length(vertex[edges[e3].v2].pos-vertex[edges[e3].v1].pos);
                    edges[e4].area = glm::length(vertex[edges[e4].v2].pos-vertex[edges[e4].v1].pos);
                    faces[f].area = resX*resY;
                }

                vertex[v1].edges[0] = e1;
                vertex[v1].edges[1] = e2;

                vertex[v2].edges[1] = e4;
                vertex[v2].edges[3] = e1;

                vertex[v3].edges[2] = e2;
                vertex[v3].edges[0] = e3;

                vertex[v4].edges[3] = e3;
                vertex[v4].edges[2] = e4;
            }
            //First Row
            else if(f/(width/resX)==0 && f%(width/resX)!=0)
            {
                edges[e1].v1 = v2;
                edges[e1].v2 = v1;
                edges[e3].v1 = v3;
                edges[e3].v2 = v4;
                edges[e4].v1 = v4;
                edges[e4].v2 = v2;

                faces[f] = Face2D(e1,-e2,e3,e4);
                if(checkVoxel(start))
                {
                    edges[e1].area = glm::length(vertex[edges[e1].v2].pos-vertex[edges[e1].v1].pos);
                    edges[e3].area = glm::length(vertex[edges[e3].v2].pos-vertex[edges[e3].v1].pos);
                    edges[e4].area = glm::length(vertex[edges[e4].v2].pos-vertex[edges[e4].v1].pos);
                    faces[f].area = resX*resY;
                }

                vertex[v2].edges[1] = e4;
                vertex[v2].edges[3] = e1;

                vertex[v3].edges[2] = e2;
                vertex[v3].edges[0] = e3;
            }
            //Other Rows First Column
            else if(f/(width/resX)!=0 && f%(width/resX)==0)
            {

                edges[e2].v1 = v1;
                edges[e2].v2 = v3;

                edges[e3].v1 = v3;
                edges[e3].v2 = v4;

                edges[e4].v1 = v4;
                edges[e4].v2 = v2;

                faces[f] = Face2D(-e1,e2,e3,e4);
                if(checkVoxel(start))
                {
                    if(y<=height/resY-1)
                    {
                        edges[e2].area = glm::length(vertex[edges[e2].v2].pos-vertex[edges[e2].v1].pos);
                    }
                    else
                    {
                        edges[e2].area = 0.0f;
                    }
                    if(y<=height/resY-1)
                    {
                        edges[e3].area = glm::length(vertex[edges[e3].v2].pos-vertex[edges[e3].v1].pos);
                    }
                    else
                    {
                        edges[e3].area = 0.0f;
                    }
                    if(y<=height/resY-1)
                    {
                        edges[e4].area = glm::length(vertex[edges[e4].v2].pos-vertex[edges[e4].v1].pos);
                    }
                    else
                    {
                        edges[e4].area = 0.0f;
                    }
                    faces[f].area = resX*resY;
                }

                vertex[v3].edges[2] = e2;
                vertex[v3].edges[0] = e3;

                vertex[v4].edges[3] = e3;
                vertex[v4].edges[2] = e4;
            }
            //Other Rows Other Column
            else if(f/(width/resX)!=0 && f%(width/resX)!=0)
            {
                edges[e3].v1 = v3;
                edges[e3].v2 = v4;

                edges[e4].v1 = v4;
                edges[e4].v2 = v2;

                faces[f] = Face2D(-e1,-e2,e3,e4);
                if(checkVoxel(start))
                {
                    if(y<=height/resY-1)
                    {
                        edges[e3].area = glm::length(vertex[edges[e3].v2].pos-vertex[edges[e3].v1].pos);
                    }
                    else
                    {
                        edges[e3].area = 0.0f;
                    }
                    if(y<=height/resY-1)
                    {
                        edges[e4].area = glm::length(vertex[edges[e4].v2].pos-vertex[edges[e4].v1].pos);
                    }
                    else
                    {
                        edges[e4].area = 0.0f;
                    }
                    faces[f].area = resX*resY;
                }

                vertex[v4].edges[3] = e3;
                vertex[v4].edges[2] = e4;
            }
            else
            {
                std::cout<<"ERROR"<<std::endl;
            }

            //0 Forward
            //1 Down
            //2 Up
            //3 Backward

/*            vertex[v1].edges[0] = e1;
            vertex[v1].edges[1] = e2;

            vertex[v2].edges[1] = e4;
            vertex[v2].edges[3] = e1;

            vertex[v3].edges[2] = e2;
            vertex[v3].edges[0] = e3;

            vertex[v4].edges[3] = e3;
            vertex[v4].edges[2] = e4;*/

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
    //derivative1(*this)*
    //Calculate Laplace Operator
    //mat = -1.0f*derivative1(*this)*hodge1(*this,true)*derivative1(*this).transpose()*hodge2(*this,false);
    std::cout<<"HODGE DIMS"<<hodge1(*this,false).cols()<<" "<<hodge1(*this,false).rows()<<std::endl;
    std::cout<<"DERIVATIVE DIMS"<<derivative1(*this,false).cols()<<" "<<derivative0(*this,false).rows()<<std::endl;

    mat = -1.0f*derivative0(*this)*hodge2(*this,true)*derivative1(*this,true)*hodge1(*this,false);
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

//TODO CALCULATE AREA OF DUAL FACES!
void Mesh2D::buildDual()
{

    dual->vertex.resize(faces.size());
    dual->edges.resize(edges.size());
    dual->faces.resize(vertex.size());
    dual->dual = this;

    for(unsigned int i=0;i<faces.size();i++)
    {
        Face2D f = faces[i];
        Edge2D e1 = edges[std::abs(f.e1)];
        Edge2D e2 = edges[std::abs(f.e2)];
        Edge2D e3 = edges[std::abs(f.e3)];
        Edge2D e4 = edges[std::abs(f.e4)];

        Vertex2D v1 = vertex[e1.v1];
        Vertex2D v2 = vertex[e1.v2];
        Vertex2D v3 = vertex[e2.v1];
        Vertex2D v4 = vertex[e2.v2];

        dual->vertex[i].pos = (1.0f/4.0f)*(v1.pos+v2.pos+v3.pos+v4.pos);
        dual->vertex[i].area = 1.0f;

        unsigned int offset = (width+resX)/resX;
        unsigned int ei1 = std::abs(f.e1);//((i*2)%(((width+resX)/resX+1)))+((i*2)/offset)*2; // Top Edge
        unsigned int ei2 = std::abs(f.e2);//((i*2)%(((width+resX)/resX+1)))+((i*2)/offset)*2+1; // Left Edge
        unsigned int ei3 = std::abs(f.e3);//((i*2)%(((width+resX)/resX+1)))+((i*2)/offset)*2+3; //Right Edge
        unsigned int ei4 = std::abs(f.e4);//((i*2)%(((width+resX)/resX+1)))+(((i+offset)*2)/offset)*2; //Bottom Edge
        //UpEdge
        if(edges[ei1].faces[0]!=-1 && edges[ei1].faces[1]!=-1 &&
           faces[edges[ei1].faces[0]].area!=0.0f && faces[edges[ei1].faces[1]].area!=0.0f)
        {
            dual->edges[ei1] = Edge2D(e1.faces[0],e1.faces[1]);
            dual->edges[ei1].area = edges[ei1].area;
        }
        else if(edges[ei1].faces[0]!=-1 && edges[ei1].faces[1]!=-1 &&
           faces[edges[ei1].faces[0]].area==0.0f && faces[edges[ei1].faces[1]].area!=0.0f)
        {
            dual->edges[ei1] = Edge2D(e1.faces[0],e1.faces[1]);
            dual->edges[ei1].area = edges[ei1].area/2.0f;
        }
        else if(edges[ei1].faces[0]!=-1 && edges[ei1].faces[1]!=-1 &&
           faces[edges[ei1].faces[0]].area!=0.0f && faces[edges[ei1].faces[1]].area==0.0f)
        {
            dual->edges[ei1] = Edge2D(e1.faces[0],e1.faces[1]);
            dual->edges[ei1].area = edges[ei1].area/2.0f;
        }
        else if((edges[ei1].faces[0]==-1 && edges[ei1].faces[1]!=-1) &&
                (faces[edges[ei1].faces[1]].area!=0.0f))
        {
            dual->edges[ei1] = Edge2D(e1.faces[0],e1.faces[1]);
            dual->edges[ei1].area = edges[ei1].area/2.0f;
        }
        else if((edges[ei1].faces[0]!=-1 && edges[ei1].faces[1]==-1) &&
                (faces[edges[ei1].faces[0]].area!=0.0f))
        {
            dual->edges[ei1] = Edge2D(e1.faces[0],e1.faces[1]);
            dual->edges[ei1].area = edges[ei1].area/2.0f;
        }
        //LeftEdge
        if(edges[ei2].faces[0]!=-1 && edges[ei2].faces[1]!=-1 &&
           faces[edges[ei2].faces[0]].area!=0.0f && faces[edges[ei2].faces[1]].area!=0.0f)
        {
            dual->edges[ei2] = Edge2D(e2.faces[0],e2.faces[1]);
            dual->edges[ei2].area = edges[ei2].area;
        }
        else if(edges[ei2].faces[0]!=-1 && edges[ei2].faces[1]!=-1 &&
           faces[edges[ei2].faces[0]].area==0.0f && faces[edges[ei2].faces[1]].area!=0.0f)
        {
            dual->edges[ei2] = Edge2D(e2.faces[0],e2.faces[1]);
            dual->edges[ei2].area = edges[ei2].area/2.0f;
        }
        else if(edges[ei2].faces[0]!=-1 && edges[ei2].faces[1]!=-1 &&
           faces[edges[ei2].faces[0]].area!=0.0f && faces[edges[ei2].faces[1]].area==0.0f)
        {
            dual->edges[ei2] = Edge2D(e2.faces[0],e2.faces[1]);
            dual->edges[ei2].area = edges[ei2].area/2.0f;
        }
        else if((edges[ei2].faces[0]==-1 && edges[ei2].faces[1]!=-1) &&
                (faces[edges[ei2].faces[1]].area!=0.0f))
        {
            dual->edges[ei2] = Edge2D(e2.faces[0],e2.faces[1]);
            dual->edges[ei2].area = edges[ei2].area/2.0f;
        }
        else if((edges[ei2].faces[0]!=-1 && edges[ei2].faces[1]==-1) &&
                (faces[edges[ei2].faces[0]].area!=0.0f))
        {
            dual->edges[ei2] = Edge2D(e2.faces[0],e2.faces[1]);
            dual->edges[ei2].area = edges[ei2].area/2.0f;
        }
        //RightEdge
        if(edges[ei3].faces[0]!=-1 && edges[ei3].faces[1]!=-1 &&
           faces[edges[ei3].faces[0]].area!=0.0f && faces[edges[ei3].faces[1]].area!=0.0f)
        {
            dual->edges[ei3] = Edge2D(e3.faces[0],e3.faces[1]);
            dual->edges[ei3].area = edges[ei3].area;
        }
        else if(edges[ei3].faces[0]!=-1 && edges[ei3].faces[1]!=-1 &&
           faces[edges[ei3].faces[0]].area==0.0f && faces[edges[ei3].faces[1]].area!=0.0f)
        {
            dual->edges[ei3] = Edge2D(e3.faces[0],e3.faces[1]);
            dual->edges[ei3].area = edges[ei3].area/2.0f;
        }
        else if(edges[ei3].faces[0]!=-1 && edges[ei3].faces[1]!=-1 &&
           faces[edges[ei3].faces[0]].area!=0.0f && faces[edges[ei3].faces[1]].area==0.0f)
        {
            dual->edges[ei3] = Edge2D(e3.faces[0],e3.faces[1]);
            dual->edges[ei3].area = edges[ei3].area/2.0f;
        }
        else if((edges[ei3].faces[0]==-1 && edges[ei3].faces[1]!=-1) &&
                (faces[edges[ei3].faces[1]].area!=0.0f))
        {
            dual->edges[ei3] = Edge2D(e3.faces[0],e3.faces[1]);
            dual->edges[ei3].area = edges[ei3].area/2.0f;
        }
        else if((edges[ei3].faces[0]!=-1 && edges[ei3].faces[1]==-1) &&
                (faces[edges[ei3].faces[0]].area!=0.0f))
        {
            dual->edges[ei3] = Edge2D(e3.faces[0],e3.faces[1]);
            dual->edges[ei3].area = edges[ei3].area/2.0f;
        }
        //DownEdge
        if(edges[ei4].faces[0]!=-1 && edges[ei4].faces[1]!=-1 &&
           faces[edges[ei4].faces[0]].area!=0.0f && faces[edges[ei4].faces[1]].area!=0.0f)
        {
            dual->edges[ei4] = Edge2D(e4.faces[0],e4.faces[1]);
            dual->edges[ei4].area = edges[ei4].area;
        }
        else if(edges[ei4].faces[0]!=-1 && edges[ei4].faces[1]!=-1 &&
           faces[edges[ei4].faces[0]].area==0.0f && faces[edges[ei4].faces[1]].area!=0.0f)
        {
            dual->edges[ei4] = Edge2D(e4.faces[0],e4.faces[1]);
            dual->edges[ei4].area = edges[ei4].area/2.0f;
        }
        else if(edges[ei4].faces[0]!=-1 && edges[ei4].faces[1]!=-1 &&
           faces[edges[ei4].faces[0]].area!=0.0f && faces[edges[ei4].faces[1]].area==0.0f)
        {
            dual->edges[ei4] = Edge2D(e4.faces[0],e4.faces[1]);
            dual->edges[ei4].area = edges[ei4].area/2.0f;
        }
        else if((edges[ei4].faces[0]==-1 && edges[ei4].faces[1]!=-1) &&
                (faces[edges[ei4].faces[1]].area!=0.0f))
        {
            dual->edges[ei4] = Edge2D(e4.faces[0],e4.faces[1]);
            dual->edges[ei4].area = edges[ei4].area/2.0f;
        }
        else if((edges[ei4].faces[0]!=-1 && edges[ei4].faces[1]==-1) &&
                (faces[edges[ei4].faces[0]].area!=0.0f))
        {
            dual->edges[ei4] = Edge2D(e4.faces[0],e4.faces[1]);
            dual->edges[ei4].area = edges[ei4].area/2.0f;
        }

        //Faces
        if(dual->edges[ei1].v1==-1&&dual->edges[ei1].v2==-1&&
           dual->edges[ei2].v1==-1&&dual->edges[ei2].v2==-1&&
           dual->edges[ei3].v1==-1&&dual->edges[ei3].v2==-1&&
           dual->edges[ei4].v1==-1&&dual->edges[ei4].v2==-1)
        {
            if(f.e1>=0)
            {
                dual->faces[e1.v2] = f;
            }
            else
            {
                dual->faces[e1.v1] = f;
            }
            std::cout<<"0";
        }
        else if(dual->edges[ei1].v2==-1)
        {
            if(dual->edges[ei4].v2==-1)
            {
                if(f.e1>=0)
                {
                    dual->faces[e1.v2] = f;
                }
                else
                {
                    dual->faces[e1.v1] = f;
                }
                std::cout<<"4";
            }
            else
            {
                if(f.e1>=0)
                {
                    dual->faces[e1.v2] = f;
                }
                else
                {
                    dual->faces[e1.v1] = f;
                }
                std::cout<<"2";
            }
        }
        else if(dual->edges[ei3].v1==-1)
        {
            if(dual->edges[ei4].v2==-1)
            {
                if(f.e1>=0)
                {
                    dual->faces[e1.v2] = f;
                }
                else
                {
                    dual->faces[e1.v1] = f;
                }
                std::cout<<"4";
            }
            else
            {
                std::cout<<"2";
            }
        }
        else if(dual->edges[ei4].v2==-1)
        {
            if(f.e1>=0)
            {
                dual->faces[e1.v2] = f;
            }
            else
            {
                dual->faces[e1.v1] = f;
            }
            std::cout<<"2";
        }
        else if(dual->edges[ei2].v2==-1)
        {
            if(f.e1>=0)
            {
                dual->faces[e1.v2] = f;
            }
            else
            {
                dual->faces[e1.v1] = f;
            }
            std::cout<<"2";
        }
        else
        {
            std::cout<<"1";
        }
        //Normal Case, add Just dual face
        if(dual->edges[ei1].v1!=-1&&dual->edges[ei1].v2!=-1&&
           dual->edges[ei2].v1!=-1&&dual->edges[ei2].v2!=-1&&
           dual->edges[ei3].v1!=-1&&dual->edges[ei3].v2!=-1&&
           dual->edges[ei4].v1!=-1&&dual->edges[ei4].v2!=-1)
        {
        }
/*
        //Up Neighbor
        if(e1.faces[0]!=-1 && e1.faces[0]!=i)
        {
        }
        if(e1.faces[1]!=-1 && e1.faces[1]!=i)
        {
        }

        //Left Neighbor
        if(e2.faces[0]!=-1 && e2.faces[0]!=i)
        {
            edges[i-2] = Edge2D(i,e4.faces[0]);
        }
        if(e2.faces[1]!=-1 && e2.faces[1]!=i)
        {
            edges[i-2] = Edge2D(i,e4.faces[1]);
        }

        //Down Neighbor
        if(e3.faces[0]!=-1 && e3.faces[0]!=i)
        {
            edges[i+1] = Edge2D(i,e4.faces[0]);
        }
        if(e3.faces[1]!=-1 && e3.faces[1]!=i)
        {
            edges[i+1] = Edge2D(i,e4.faces[1]);
        }

        //Right Neighbor
        if(e4.faces[0]!=-1 && e4.faces[0]!=i)
        {
            edges[i] = Edge2D(i,e4.faces[1]);
        }
        if(e4.faces[1]!=-1 && e4.faces[1]!=i)
        {
            edges[i] = Edge2D(i,e4.faces[1]);
        }
        */
    }
}

