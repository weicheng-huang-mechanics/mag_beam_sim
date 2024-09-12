#ifndef EXTERNALMAGNETICFORCE_H
#define EXTERNALMAGNETICFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalMagneticForce
{
public:
	externalMagneticForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~externalMagneticForce();

	void computeFm(double m_time);
	void computeJm();

private:
	elasticPlate *plate;
	timeStepper *stepper;

    int currentIndex;

    Vector2d baVector;
    Vector2d brVector;

    Matrix2d baGradient;

    Vector2d m_current;
    Vector2d t_current;

    Vector2d m_start;
    Vector2d t_start;

    Matrix2d Id2;

    Matrix2d dtde;
    Matrix2d dmde;
    Vector2d dEde;

    Vector2d x1;
    Vector2d x2;

    Vector2d x1_start;
    Vector2d x2_start;

    double edge;

    VectorXd force;
    VectorXi localDOF;

    Vector2d M;

    double rotateAngle;

    int nv1;
    int nv2;
};

#endif
