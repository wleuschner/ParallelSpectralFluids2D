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
    for(unsigned int i=0;i<mesh->getNumFaces();i++)
    {
    }
}
