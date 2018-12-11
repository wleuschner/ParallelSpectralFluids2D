#define GLM_ENABLE_EXPERIMENTAL

#include "flowviswidget.h"
#include <glm/glm.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <Eigen/Eigen>
#include <QPainter>


FlowVisWidget::FlowVisWidget(QWidget *parent) : QWidget(parent)
{

}

void FlowVisWidget::setMesh(Mesh2D* mesh)
{
    this->mesh = mesh;
}

void FlowVisWidget::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    bool* voxelMap = mesh->getVoxelMap();
    unsigned int width = mesh->getWidth();
    unsigned int height = mesh->getHeight();
    unsigned int resX = mesh->getResX();
    unsigned int resY = mesh->getResY();

    for(std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>>::iterator fit=mesh->faces.begin();fit!=mesh->faces.end();fit++)
    {
        QPen pen(Qt::black, 2, Qt::SolidLine);
        painter.setPen(pen);
        unsigned int iv1 = std::get<0>(*fit);
        unsigned int iv2 = std::get<1>(*fit);
        unsigned int iv3 = std::get<2>(*fit);
        unsigned int iv4 = std::get<3>(*fit);


        Vertex2D v1 = mesh->vertex[iv1];
        Vertex2D v2 = mesh->vertex[iv2];
        Vertex2D v3 = mesh->vertex[iv3];
        Vertex2D v4 = mesh->vertex[iv4];

        unsigned int s1=1;
        unsigned int s2=1;
        unsigned int s3=1;
        unsigned int s4=1;
        std::tuple<unsigned int,unsigned int> e1,e2,e3,e4;
        std::set<std::tuple<unsigned int,unsigned int>>::iterator eit1,eit2,eit3,eit4;
        if((eit1=mesh->edges.find(std::make_tuple(iv1,iv2)))==mesh->edges.end())
        {
            s1=-1;
            eit1 = mesh->edges.find(std::make_tuple(iv2,iv1));
        }
        if((eit2=mesh->edges.find(std::make_tuple(iv2,iv3)))==mesh->edges.end())
        {
            s2=-1;
            eit2 = mesh->edges.find(std::make_tuple(iv3,iv2));
        }
        if((eit3=mesh->edges.find(std::make_tuple(iv3,iv4)))==mesh->edges.end())
        {
            s3=-1;
            eit3 = mesh->edges.find(std::make_tuple(iv4,iv3));
        }
        if((eit4=mesh->edges.find(std::make_tuple(iv4,iv1)))==mesh->edges.end())
        {
            s4=-1;
            eit4 = mesh->edges.find(std::make_tuple(iv1,iv4));
        }
        glm::vec2 vel = glm::rotate(s1*mesh->velocityField(std::distance(mesh->edges.begin(),eit1))*(v2.pos-v1.pos),90.0f)+
                        glm::rotate(s2*mesh->velocityField(std::distance(mesh->edges.begin(),eit2))*(v3.pos-v2.pos),90.0f)+
                        glm::rotate(s3*mesh->velocityField(std::distance(mesh->edges.begin(),eit3))*(v4.pos-v3.pos),90.0f)+
                        glm::rotate(s4*mesh->velocityField(std::distance(mesh->edges.begin(),eit4))*(v1.pos-v4.pos),90.0f);
        vel = 0.000000000002f*vel;

        painter.drawLine(v1.pos.x,v1.pos.y,v2.pos.x,v2.pos.y);
        painter.drawLine(v2.pos.x,v2.pos.y,v3.pos.x,v3.pos.y);
        painter.drawLine(v3.pos.x,v3.pos.y,v4.pos.x,v4.pos.y);
        painter.drawLine(v4.pos.x,v4.pos.y,v1.pos.x,v1.pos.y);

        painter.drawLine(v3.pos.x+(v1.pos.x-v3.pos.x)/2,v3.pos.y+(v1.pos.y-v3.pos.y)/2,v3.pos.x+(v1.pos.x-v3.pos.x)/2+vel.x,v3.pos.y+(v1.pos.y-v3.pos.y)/2+vel.y);
        painter.drawEllipse(v3.pos.x+(v1.pos.x-v3.pos.x)/2-4,v3.pos.y+(v1.pos.y-v3.pos.y)/2-4,8.0,8.0);
    }
    /*
    for(unsigned int y=0;y<height/resY;y++)
    {
        for(unsigned int x=0;x<width/resX;x++)
        {
            if(voxelMap[y*(width/resX)+x])
            {
                QPen pen(Qt::black, 2, Qt::SolidLine);
                painter.setPen(pen);
                painter.drawRect(x*resX,y*resY,resX,resY);
                Face2D face = mesh->faces[y*(width/resX)+x];
                glm::vec2 dir(0.0f,0.0f);
                for(unsigned int e=0;e<4;e++)
                {
                    Edge2D edge = mesh->edges[std::abs(face.e[e])];
                    Vertex2D v1 = mesh->vertex[edge.v1];
                    Vertex2D v2 = mesh->vertex[edge.v2];
                    if(face.e[e]>=0)
                    {
                        dir+=mesh->velocityField(std::abs(face.e[e]))*(v1.pos.x-v2.pos);
                    }
                    else
                    {
                        dir+=mesh->velocityField(std::abs(face.e[e]))*(v2.pos.x-v1.pos);
                    }
                }
                pen.setColor(QColor(255*glm::length(dir),0,0));
                painter.setPen(pen);
                painter.drawLine(x*resX+0.5f*resX,y*resY+0.5f*resY,x*resX+0.5f*resX+dir.x*100,y*resY+0.5f*resY+dir.y*100);
            }
        }
    }*/
}
