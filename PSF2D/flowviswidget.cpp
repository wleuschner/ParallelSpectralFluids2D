#include "flowviswidget.h"
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

FlowVisWidget::FlowVisWidget(QWidget *parent) : QOpenGLWidget(parent)
{
    velocityNormalizationState = UNIT_NORMALIZATION;
    meshAlpha = 0.3;
    gridVisible = true;
    velocityVisible = true;
    vorticityVisible = false;
    recording = false;
    imageNo = 0;
}

void FlowVisWidget::setSolver(AbstractSolver* solver)
{
    this->solver = solver;
    resizeGrid();
}

void FlowVisWidget::resizeGrid()
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

        glm::vec2 e1 = v2.pos-v1.pos;
        glm::vec2 e2 = v3.pos-v2.pos;
        glm::vec2 e3 = v4.pos-v3.pos;
        glm::vec2 e4 = v1.pos-v4.pos;

        glm::vec2 center = v3.pos+0.5f*(v1.pos-v3.pos);

        gridEdges[ie1] = e1;
        gridEdges[ie2] = e2;
        gridEdges[ie3] = e3;
        gridEdges[ie4] = e4;

        faceCenters[decMesh.getFaceIndex(iv1,iv2,iv3,iv4)] = center;
    }
}

void FlowVisWidget::setImage(QImage image)
{
    originalImage = image;
    rotationField = QImage(originalImage.width(),originalImage.height(),QImage::Format_RGBA8888);
    renderBuffer = QImage(originalImage.width(),originalImage.height(),QImage::Format_RGBA8888);
    rotationFieldPixels = (unsigned int*)rotationField.bits();

}

void FlowVisWidget::initializeGL()
{
    program.addShaderFromSourceFile(QOpenGLShader::Vertex,"../PSF2D/flow.vs");
    program.addShaderFromSourceFile(QOpenGLShader::Fragment,"../PSF2D/flow.fs");
    program.link();

    glClearColor(1.0f,1.0f,1.0f,1.0f);
    glDisable(GL_DEPTH_TEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
    glShadeModel(GL_SMOOTH);
}

void FlowVisWidget::resizeGL(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,width(),height(),0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void FlowVisWidget::setVelocityNormalization(VelocityNormalization state)
{
    velocityNormalizationState = state;
}

void FlowVisWidget::setGridVisisble(bool state)
{
    gridVisible = state;
}

void FlowVisWidget::setVelocityVisible(bool state)
{
    velocityVisible = state;
}

void FlowVisWidget::setVorticityVisible(bool state)
{
    vorticityVisible = state;
}


/*void FlowVisWidget::paintGL()
{
    DECMesh2D decMesh = solver->getDECMesh();
    Eigen::VectorXf velocityField = solver->getVelocityField();
    Eigen::VectorXf vorticityField = solver->getVorticityField();
    unsigned int width = solver->getMesh()->getWidth();
    unsigned int height = solver->getMesh()->getHeight();
    unsigned int resX = solver->getMesh()->getResolution();
    unsigned int resY = solver->getMesh()->getResolution();

    float maxVort = solver->getMaxVorticity();
//    float minVort = mesh->vorticityField.cwiseAbs().minCoeff();
    //float maxVort = mesh->maxRotation;
    float minVort = solver->getMinVorticity();

    for(FaceIterator fit=solver->getDECMesh().getFaceIteratorBegin();fit!=solver->getDECMesh().getFaceIteratorEnd();fit++)
    {

        unsigned int iv1 = std::get<0>(*fit);
        unsigned int iv2 = std::get<1>(*fit);
        unsigned int iv3 = std::get<2>(*fit);
        unsigned int iv4 = std::get<3>(*fit);


        Vertex2D v1 = solver->getMesh()->vertex[iv1];
        Vertex2D v2 = solver->getMesh()->vertex[iv2];
        Vertex2D v3 = solver->getMesh()->vertex[iv3];
        Vertex2D v4 = solver->getMesh()->vertex[iv4];

        int s1=decMesh.getEdgeSignum(iv1,iv2);
        int s2=decMesh.getEdgeSignum(iv2,iv3);
        int s3=decMesh.getEdgeSignum(iv3,iv4);
        int s4=decMesh.getEdgeSignum(iv4,iv1);

        if(vorticityVisible)
        {
            glm::vec4 vorts((vorticityField(decMesh.getPointIndex(iv1)))/(maxVort*0.5),
                            (vorticityField(decMesh.getPointIndex(iv2)))/(maxVort*0.5),
                            (vorticityField(decMesh.getPointIndex(iv3)))/(maxVort*0.5),
                            (vorticityField(decMesh.getPointIndex(iv4)))/(maxVort*0.5));
            vorts = glm::clamp(vorts,glm::vec4(-1.0f,-1.0f,-1.0f,-1.0f),glm::vec4(1.0f,1.0f,1.0f,1.0f));
            vorts += glm::vec4(1.0f,1.0f,1.0f,1.0f);

            glm::vec3 color1 = glm::mix(glm::vec3(1.0f,1.0f,1.0f),glm::vec3(0.0f,0.0f,0.0f),vorts.x/2);
            glm::vec3 color2 = glm::mix(glm::vec3(1.0f,1.0f,1.0f),glm::vec3(0.0f,0.0f,0.0f),vorts.y/2);
            glm::vec3 color3 = glm::mix(glm::vec3(1.0f,1.0f,1.0f),glm::vec3(0.0f,0.0f,0.0f),vorts.z/2);
            glm::vec3 color4 = glm::mix(glm::vec3(1.0f,1.0f,1.0f),glm::vec3(0.0f,0.0f,0.0f),vorts.w/2);


            glLineWidth(1.0f);
            glBegin(GL_LINES);
            for(unsigned int y=0;y<resY;y++)
            {
                float yOfs = static_cast<float>(y)/static_cast<float>(resY);
                glm::vec2 a=glm::mix(glm::vec2(vort1,vort4),glm::vec2(vort2,vort3),yOfs);
                a = glm::clamp(a,glm::vec2(-1.0f,-1.0f),glm::vec2(1.0f,1.0f));
                a += glm::vec2(1.0f,1.0f);
                glm::vec3 color1 = glm::mix(glm::vec3(1.0f,1.0f,1.0f),glm::vec3(0.0f,0.0f,0.0f),a.x/2.0f);
                glm::vec3 color2 = glm::mix(glm::vec3(1.0f,1.0f,1.0f),glm::vec3(0.0f,0.0f,0.0f),a.y/2.0f);
                glColor3ub(color1.r*255,color1.g*255,color1.b*255);
                glVertex2f(v1.pos.x,v1.pos.y+y);
                glColor3ub(color2.r*255,color2.g*255,color2.b*255);
                glVertex2f(v4.pos.x,v4.pos.y+y);
            }
            glEnd();
            glLineWidth(2.0f);


            program.bind();
            glBegin(GL_QUADS);
            glColor3f(color1.r,color1.g,color1.b);
            glVertex2f(v1.pos.x,v1.pos.y);
            glColor3f(color2.r,color2.g,color2.b);
            glVertex2f(v2.pos.x,v2.pos.y);
            glColor3f(color3.r,color3.g,color3.b);
            glVertex2f(v3.pos.x,v3.pos.y);
            glColor3f(color4.r,color4.g,color4.b);
            glVertex2f(v4.pos.x,v4.pos.y);
            glEnd();
            program.release();
        }

        if(gridVisible)
        {
            glBegin(GL_LINES);
            glColor3f(0.0f,0.0f,0.0f);
            glVertex2f(v1.pos.x,v1.pos.y);
            glVertex2f(v2.pos.x,v2.pos.y);
            glVertex2f(v2.pos.x,v2.pos.y);
            glVertex2f(v3.pos.x,v3.pos.y);
            glVertex2f(v3.pos.x,v3.pos.y);
            glVertex2f(v4.pos.x,v4.pos.y);
            glVertex2f(v4.pos.x,v4.pos.y);
            glVertex2f(v1.pos.x,v1.pos.y);
            glEnd();
        }

        if(velocityVisible)
        {
            glm::vec2 e1 = v2.pos-v1.pos;
            glm::vec2 e2 = v3.pos-v2.pos;
            glm::vec2 e3 = v4.pos-v3.pos;
            glm::vec2 e4 = v1.pos-v4.pos;

            glm::vec2 vel = glm::rotate(s1*velocityField(decMesh.getEdgeIndex(iv1,iv2))*glm::normalize(e1)+
                                        s2*velocityField(decMesh.getEdgeIndex(iv2,iv3))*glm::normalize(e2)+
                                        s3*velocityField(decMesh.getEdgeIndex(iv3,iv4))*glm::normalize(e3)+
                                        s4*velocityField(decMesh.getEdgeIndex(iv4,iv1))*glm::normalize(e4),glm::radians(-90.0f));
            vel = solver->getTimestep()*vel;

            glm::vec2 center = v3.pos+0.5f*(v1.pos-v3.pos);
            glm::vec2 dir = center+vel;
            glBegin(GL_LINES);
            glColor3f(1.0f,0.0f,0.0f);
            glVertex2f(center.x,center.y);
            glVertex2f(dir.x,dir.y);
            glEnd();
        }
    }
    if(recording)
    {

        grabFramebuffer().save("image"+QString::number(imageNo)+".bmp","bmp");
        imageNo++;
    }
}*/

void FlowVisWidget::paintEvent(QPaintEvent *event)
{

    //QImage rotationField(originalImage.width(),originalImage.height(),QImage::Format_RGBA8888);


    QPainter painter;
    DECMesh2D decMesh = solver->getDECMesh();
    Eigen::VectorXf velocityField = solver->getVelocityField();
    Eigen::VectorXf vorticityField = solver->getVorticityField();
    unsigned int width = solver->getMesh()->getWidth();
    unsigned int height = solver->getMesh()->getHeight();
    unsigned int resX = solver->getMesh()->getResolution();
    unsigned int resY = solver->getMesh()->getResolution();

    renderBuffer.fill(Qt::transparent);
    float maxVort = solver->getMaxVorticity();
//    float minVort = mesh->vorticityField.cwiseAbs().minCoeff();
    //float maxVort = mesh->maxRotation;
    float minVort = solver->getMinVorticity();

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

                float vort1 = (vorticityField(decMesh.getPointIndex(iv1)))/(maxVort*0.05);
                float vort2 = (vorticityField(decMesh.getPointIndex(iv2)))/(maxVort*0.05);
                float vort3 = (vorticityField(decMesh.getPointIndex(iv3)))/(maxVort*0.05);
                float vort4 = (vorticityField(decMesh.getPointIndex(iv4)))/(maxVort*0.05);

                glm::vec2 ySlope = (1.0f/resY)*glm::vec2(vort2-vort1,vort3-vort4);
                glm::vec2 yBegin = glm::vec2(vort1,vort4);
                for(unsigned int y=0;y<resY;y++)
                {
                    unsigned int posY = v1.pos.y+y;
                    unsigned int *startLine = rotationFieldPixels+((posY)*imageWidth+static_cast<unsigned int>(v1.pos.x));

                    float xSlope = (1.0f/resX)*(yBegin.y-yBegin.x);
                    float xBegin = yBegin.x;
                    for(unsigned int x=0;x<resX;x++)
                    {
                        float c = xBegin;
                        c = glm::clamp(c,-1.0f,1.0f);
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

                glm::vec2 e1 = gridEdges[ie1];
                glm::vec2 e2 = gridEdges[ie2];
                glm::vec2 e3 = gridEdges[ie3];
                glm::vec2 e4 = gridEdges[ie4];

                float s1=decMesh.getEdgeSignum(iv1,iv2);
                float s2=decMesh.getEdgeSignum(iv2,iv3);
                float s3=decMesh.getEdgeSignum(iv3,iv4);
                float s4=decMesh.getEdgeSignum(iv4,iv1);

                glm::vec2 vel = glm::rotate(0.5f*(s1*velocityField(ie1)*glm::normalize(s1*e1)+s3*velocityField(ie3)*glm::normalize(s3*e3))+
                                            0.5f*(s2*velocityField(ie2)*glm::normalize(s2*e2)+s4*velocityField(ie4)*glm::normalize(s4*e4)),
                                            glm::radians(-90.0f));
                if(velocityNormalizationState==UNIT_NORMALIZATION)
                {
                    vel = (solver->getMesh()->getResolution()/2.0f)*glm::normalize(vel);
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
            glm::vec2 faceCenter = faceCenters[fidx];
            glm::vec2 vel = velocities[fidx];

            painter.drawEllipse(faceCenter.x-4,faceCenter.y-4,8.0,8.0);

            glm::vec2 dir = faceCenter+vel;
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

void FlowVisWidget::keyPressEvent(QKeyEvent *event)
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
