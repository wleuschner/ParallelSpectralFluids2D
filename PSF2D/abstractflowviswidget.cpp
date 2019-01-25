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
    DECMesh2D decMesh = solver->getDECMesh();
    velocities.resize(decMesh.getNumFaces());
    faceCenters.resize(decMesh.getNumFaces());
    gridEdges.resize(decMesh.getNumEdges());


    for(FaceIterator fit=solver->getDECMesh().getFaceIteratorBegin();fit!=solver->getDECMesh().getFaceIteratorEnd();fit++)
    {
        unsigned int iv1 = std::get<0>(*fit);
        unsigned int iv2 = std::get<1>(*fit);
        unsigned int iv3 = std::get<2>(*fit);
        unsigned int iv4 = std::get<3>(*fit);

        unsigned int ie1 = decMesh.getEdgeIndex(iv1,iv2);
        unsigned int ie2 = decMesh.getEdgeIndex(iv2,iv3);
        unsigned int ie3 = decMesh.getEdgeIndex(iv3,iv4);
        unsigned int ie4 = decMesh.getEdgeIndex(iv4,iv1);

        Vertex2D v1 = solver->getMesh()->vertex[iv1];
        Vertex2D v2 = solver->getMesh()->vertex[iv2];
        Vertex2D v3 = solver->getMesh()->vertex[iv3];
        Vertex2D v4 = solver->getMesh()->vertex[iv4];

        glm::dvec2 e1 = v2.pos-v1.pos;
        glm::dvec2 e2 = v3.pos-v2.pos;
        glm::dvec2 e3 = v4.pos-v3.pos;
        glm::dvec2 e4 = v1.pos-v4.pos;

        glm::dvec2 center = v3.pos+0.5*(v1.pos-v3.pos);

        gridEdges[ie1] = e1;
        gridEdges[ie2] = e2;
        gridEdges[ie3] = e3;
        gridEdges[ie4] = e4;

        faceCenters[decMesh.getFaceIndex(iv1,iv2,iv3,iv4)] = center;
    }
}

void AbstractFlowVisWidget::setImage(QImage image)
{
    originalImage = image;
    rotationField = QImage(originalImage.width(),originalImage.height(),QImage::Format_RGBA8888);
    renderBuffer = QImage(originalImage.width(),originalImage.height(),QImage::Format_RGBA8888);
    rotationFieldPixels = (unsigned int*)rotationField.bits();

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

    //QImage rotationField(originalImage.width(),originalImage.height(),QImage::Format_RGBA8888);


    QPainter painter;
    DECMesh2D decMesh = solver->getDECMesh();
    Eigen::VectorXd velocityField = solver->getVelocityField();
    Eigen::VectorXd vorticityField = solver->getVorticityField();
    unsigned int width = solver->getMesh()->getWidth();
    unsigned int height = solver->getMesh()->getHeight();
    unsigned int resX = solver->getMesh()->getResolution();
    unsigned int resY = solver->getMesh()->getResolution();

    renderBuffer.fill(Qt::transparent);
    double maxVort = solver->getMaxVorticity();
//    double minVort = mesh->vorticityField.cwiseAbs().minCoeff();
    //double maxVort = mesh->maxRotation;
    double minVort = solver->getMinVorticity();

    painter.begin(&renderBuffer);

    QPen black(QColor(0,0,0), 2, Qt::SolidLine);
    QPen red(QColor(255,0,0), 2, Qt::SolidLine);

    rotationField.fill(Qt::white);
    unsigned int imageWidth = originalImage.width();
    unsigned int imageHeight = originalImage.height();

    #pragma omp parallel
    {
        FaceIterator fit=solver->getDECMesh().getFaceIteratorBegin();
        std::advance(fit,omp_get_thread_num()*(solver->getDECMesh().getNumFaces()/omp_get_num_threads()));
        #pragma omp for
        for(unsigned int i=0;i<decMesh.getNumFaces();i++)
        {
            unsigned int iv1 = std::get<0>(*fit);
            unsigned int iv2 = std::get<1>(*fit);
            unsigned int iv3 = std::get<2>(*fit);
            unsigned int iv4 = std::get<3>(*fit);

            if(vorticityVisible)
            {
                Vertex2D v1 = solver->getMesh()->vertex[iv1];

                double vort1 = (vorticityField(decMesh.getPointIndex(iv1)))/(maxVort*0.05);
                double vort2 = (vorticityField(decMesh.getPointIndex(iv2)))/(maxVort*0.05);
                double vort3 = (vorticityField(decMesh.getPointIndex(iv3)))/(maxVort*0.05);
                double vort4 = (vorticityField(decMesh.getPointIndex(iv4)))/(maxVort*0.05);

                glm::dvec2 ySlope = (1.0/resY)*glm::dvec2(vort2-vort1,vort3-vort4);
                glm::dvec2 yBegin = glm::dvec2(vort1,vort4);
                for(unsigned int y=0;y<resY;y++)
                {
                    unsigned int posY = v1.pos.y+y;
                    unsigned int *startLine = rotationFieldPixels+((posY)*imageWidth+static_cast<unsigned int>(v1.pos.x));

                    double xSlope = (1.0f/resX)*(yBegin.y-yBegin.x);
                    double xBegin = yBegin.x;
                    for(unsigned int x=0;x<resX;x++)
                    {
                        double c = xBegin;
                        c = glm::clamp(c,-1.0,1.0);
                        glm::ivec3 color;
                        if(c>0)
                        {
                            color = glm::mix(glm::ivec3(255,255,255),glm::ivec3(0,0,255),c);
                        }
                        else
                        {
                            color = glm::mix(glm::ivec3(255,255,255),glm::ivec3(0,255,0),-c);
                        }
                        *startLine = qRgb(color.b,color.g,color.r);
                        startLine++;
                        xBegin += xSlope;
                    }
                    yBegin += ySlope;
                }
            }

            if(velocityVisible)
            {
                unsigned int ie1 = decMesh.getEdgeIndex(iv1,iv2);
                unsigned int ie2 = decMesh.getEdgeIndex(iv2,iv3);
                unsigned int ie3 = decMesh.getEdgeIndex(iv3,iv4);
                unsigned int ie4 = decMesh.getEdgeIndex(iv4,iv1);

                unsigned int fidx = decMesh.getFaceIndex(iv1,iv2,iv3,iv4);

                glm::dvec2 e1 = gridEdges[ie1];
                glm::dvec2 e2 = gridEdges[ie2];
                glm::dvec2 e3 = gridEdges[ie3];
                glm::dvec2 e4 = gridEdges[ie4];

                double s1=decMesh.getEdgeSignum(iv1,iv2);
                double s2=decMesh.getEdgeSignum(iv2,iv3);
                double s3=decMesh.getEdgeSignum(iv3,iv4);
                double s4=decMesh.getEdgeSignum(iv4,iv1);

                glm::dvec2 vel = glm::rotate(0.5*(s1*velocityField(ie1)*glm::normalize(s1*e1)+s3*velocityField(ie3)*glm::normalize(s3*e3))+
                                            0.5*(s2*velocityField(ie2)*glm::normalize(s2*e2)+s4*velocityField(ie4)*glm::normalize(s4*e4)),
                                            glm::radians(-90.0));
                if(velocityNormalizationState==UNIT_NORMALIZATION)
                {
                    vel = (solver->getMesh()->getResolution()/2.0)*glm::normalize(vel);
                }
                else if(velocityNormalizationState==TIMESTEP_NORMALIZATION)
                {
                    vel = solver->getTimestep()*vel;
                }

                velocities[fidx] = vel;
            }
            std::advance(fit,1);
        }
    }

    for(FaceIterator fit=solver->getDECMesh().getFaceIteratorBegin();fit!=solver->getDECMesh().getFaceIteratorEnd();fit++)
    {
        painter.setPen(black);
        unsigned int iv1 = std::get<0>(*fit);
        unsigned int iv2 = std::get<1>(*fit);
        unsigned int iv3 = std::get<2>(*fit);
        unsigned int iv4 = std::get<3>(*fit);

        unsigned int fidx = decMesh.getFaceIndex(iv1,iv2,iv3,iv4);

        Vertex2D v1 = solver->getMesh()->vertex[iv1];
        Vertex2D v2 = solver->getMesh()->vertex[iv2];
        Vertex2D v3 = solver->getMesh()->vertex[iv3];
        Vertex2D v4 = solver->getMesh()->vertex[iv4];

        if(gridVisible)
        {
            painter.drawLine(v1.pos.x,v1.pos.y,v2.pos.x,v2.pos.y);
            painter.drawLine(v2.pos.x,v2.pos.y,v3.pos.x,v3.pos.y);
            painter.drawLine(v3.pos.x,v3.pos.y,v4.pos.x,v4.pos.y);
            painter.drawLine(v4.pos.x,v4.pos.y,v1.pos.x,v1.pos.y);
        }

        if(velocityVisible)
        {
            glm::dvec2 faceCenter = faceCenters[fidx];
            glm::dvec2 vel = velocities[fidx];

            painter.drawEllipse(faceCenter.x-4,faceCenter.y-4,8.0,8.0);

            glm::dvec2 dir = faceCenter+vel;
            painter.setPen(red);
            painter.drawLine(faceCenter.x,faceCenter.y,dir.x,dir.y);
        }
    }
    painter.end();

    QImage image(originalImage.width(),originalImage.height(),QImage::Format_RGBA8888);
    painter.begin(&image);
    painter.drawImage(0,0,rotationField);
    painter.drawImage(0,0,renderBuffer);
    painter.end();

    painter.begin(this);
    painter.drawImage(0,0,image);
    painter.end();
    if(recording)
    {

        image.save("image"+QString::number(imageNo)+".bmp","bmp");
        imageNo++;
    }
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
