#include "voxelviswidget.h"
#include <QPainter>

VoxelVisWidget::VoxelVisWidget(QWidget *parent) : QWidget(parent)
{

}

void VoxelVisWidget::setMesh(Mesh2D* mesh)
{
    this->mesh = mesh;
}

void VoxelVisWidget::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    bool* voxelMap = mesh->getVoxelMap();
    unsigned int width = mesh->getWidth();
    unsigned int height = mesh->getHeight();
    unsigned int resX = mesh->getResX();
    unsigned int resY = mesh->getResY();
    for(unsigned int i=0;i<mesh->getNumEdges();i++)
    {
        Edge2D edge = mesh->edges[i];
        if(edge.area!=0.0f)
        {
            Vertex2D v1 = mesh->vertex[edge.v1];
            Vertex2D v2 = mesh->vertex[edge.v2];
            QPen pen(Qt::black, 2, Qt::SolidLine);
            painter.setPen(pen);
            //painter.drawPoint(v1.pos.x,v1.pos.y);
            //painter.drawPoint(v2.pos.x,v2.pos.y);
            painter.drawLine(v1.pos.x,v1.pos.y,v2.pos.x,v2.pos.y);
        }
    }
}
