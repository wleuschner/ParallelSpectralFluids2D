#ifndef VOXELVISWIDGET_H
#define VOXELVISWIDGET_H

#include <QWidget>
#include "mesh2d.h"

class VoxelVisWidget : public QWidget
{
    Q_OBJECT
public:
    explicit VoxelVisWidget(QWidget *parent = nullptr);
    void setMesh(Mesh2D* mesh);
protected:
    void paintEvent(QPaintEvent *event);
signals:

public slots:
private:
    Mesh2D* mesh;
};

#endif // VOXELVISWIDGET_H
