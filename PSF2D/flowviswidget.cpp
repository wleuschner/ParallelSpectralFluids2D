#include "flowviswidget.h"
#include <glm/glm.hpp>
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
    for(unsigned int y=0;y<height/resY;y++)
    {
        for(unsigned int x=0;x<width/resX;x++)
        {
/*            if(voxelMap[y*(width/resX)+x])
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
            }*/
        }
    }
}
