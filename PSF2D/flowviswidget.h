#ifndef FLOWVISWIDGET_H
#define FLOWVISWIDGET_H

#include <QObject>
#include <QWidget>
#include "mesh2d.h"

class FlowVisWidget : public QWidget
{
    Q_OBJECT
public:
    explicit FlowVisWidget(QWidget *parent = nullptr);
    void setMesh(Mesh2D* mesh);
protected:
    void paintEvent(QPaintEvent *event);
signals:

public slots:
private:
    Mesh2D* mesh;
};

#endif // FLOWVISWIDGET_H
