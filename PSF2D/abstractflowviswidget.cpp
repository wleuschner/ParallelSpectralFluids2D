#include "abstractflowviswidget.h"
#include <glm/glm.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <Eigen/Eigen>
#include <QKeyEvent>
#include <QPixmap>
#include <QPainter>
#include <iostream>
#include <GL/glu.h>
#include <QOpenGLShader>
#include <omp.h>

AbstractFlowVisWidget::AbstractFlowVisWidget(QWidget *parent) : QOpenGLWidget(parent)
{
    velocityNormalizationState = UNIT_NORMALIZATION;
    meshAlpha = 0.3;
    gridVisible = true;
    velocityVisible = true;
    vorticityVisible = false;
    recording = false;
    imageNo = 0;
}

void AbstractFlowVisWidget::setSolver(AbstractSolver* solver)
{
    this->solver = solver;
    resizeGrid();
}

void AbstractFlowVisWidget::resizeGrid()
{

}

void AbstractFlowVisWidget::setImage(QImage image)
{

}

void AbstractFlowVisWidget::initializeGL()
{
    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glDisable(GL_DEPTH_TEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
    glShadeModel(GL_SMOOTH);
}

void AbstractFlowVisWidget::resizeGL(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,width(),height(),0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void AbstractFlowVisWidget::setVelocityNormalization(VelocityNormalization state)
{
    velocityNormalizationState = state;
}

void AbstractFlowVisWidget::setGridVisisble(bool state)
{
    gridVisible = state;
}

void AbstractFlowVisWidget::setVelocityVisible(bool state)
{
    velocityVisible = state;
}

void AbstractFlowVisWidget::setVorticityVisible(bool state)
{
    vorticityVisible = state;
}

void AbstractFlowVisWidget::paintEvent(QPaintEvent *event)
{

}

void AbstractFlowVisWidget::keyPressEvent(QKeyEvent *event)
{
    switch(event->key())
    {
        case Qt::Key_Space:
        {
            if(recording)
            {
                recording = false;
                imageNo=0;
            }
            else
            {
                recording = true;
            }
        }
    }
}
