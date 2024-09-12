#include "externalMagneticForce.h"

externalMagneticForce::externalMagneticForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;

	baGradient.setZero(2, 2);

	Id2<<1,0,
         0,1;

    force.setZero(4, 1);
    localDOF.setZero(4, 1);

    currentIndex = 0;
}

externalMagneticForce::~externalMagneticForce()
{
	;
}

void externalMagneticForce::computeFm(double m_time)
{
	if (currentIndex < plate->timeSeries.size() - 1 && plate->timeSeries[currentIndex + 1] <= m_time)
	{
		currentIndex++;
	}

	baVector = plate->ba_vec[currentIndex];

	for (int i = 0; i < plate->edgeNum; i++)
	{
		nv1 = plate->v_edgeElement[i].nv_1;
		nv2 = plate->v_edgeElement[i].nv_2;

		brVector = plate->br_vec[i];

		x1 = plate->getVertexOld(nv1); 
		x2 = plate->getVertexOld(nv2); 


		t_current = (x2 - x1) / (x2 - x1).norm();
		m_current(0) = - t_current(1);
		m_current(1) = t_current(0);

		edge = (x2 - x1).norm();

		x1_start = plate->getVertexStart(nv1); 
		x2_start = plate->getVertexStart(nv2); 

		t_start = (x2_start - x1_start) / (x2_start - x1_start).norm();
		m_start(0) = - t_start(1);
		m_start(1) = t_start(0);

		dtde = (Id2 - t_current * t_current.transpose()) / edge;
		dmde = - (t_current * m_current.transpose()) / edge;

		M = m_start.dot(brVector) * m_current + t_start.dot(brVector) * t_current;

		dEde = (m_start.dot(brVector) * dmde + t_start.dot(brVector) * dtde).transpose() * baVector;

		force.segment(0, 2) = - dEde + baGradient * M / 2;
		force.segment(2, 2) = dEde + baGradient * M / 2;

		//cout << i << " " << force.norm() << endl;

		force = - force * edge * plate->crossSectionalArea;

		localDOF(0) = 2 * nv1 + 0;
		localDOF(1) = 2 * nv1 + 1;
		localDOF(2) = 2 * nv2 + 0;
		localDOF(3) = 2 * nv2 + 1;

		for (int k = 0; k < 4; k++)
		{
			stepper->addForce(localDOF(k), - force(k));
		}

	}
}

void externalMagneticForce::computeJm()
{
	;
}