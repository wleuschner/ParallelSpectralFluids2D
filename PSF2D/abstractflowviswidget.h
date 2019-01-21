#ifndef ABSTRACTFLOWVISWIDGET_H
#define ABSTRACTFLOWVISWIDGET_H

#include <QObject>
#include <QOpenGLWidget>

class AbstractFlowVisWidget : public QOpenGLWidget
{
    Q_OBJECT
public:
    explicit AbstractFlowVisWidget(QWidget *parent = nullptr);

signals:

public slots:
};

#endif // ABSTRACTFLOWVISWIDGET_H
