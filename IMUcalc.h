#include <math.h>
#include "BOARD.h"

double psiFromDCM(double DCMMatrix[3][3]);

double thetaFromDCM(double DCMMatrix[3][3]);

double phiFromDCM(double DCMMatrix[3][3]);

void EulFromDCM(double DCMMatrix[3][3],double Euler[3]);

double sinc(double x);

void Rexp(double w[3], double ans[3][3]);

double norm(double x[3]);

void IntegrateOpenLoopMatExp(double Rminus[3][3], double gyro[3], double deltaT, double Rplus[3][3]);

void IntegrateOpenLoopForwardInt(double Rminus[3][3], double gyro[3], double deltaT, double Rplus[3][3]);

void IntegrateClosedLoop(double Rminus[3][3], double Bminus[3], double gyros[3], double mags[3], double accel[3], double dT, double Rplus[3][3], double Bplus[3]);

void DCMfromTriad(double accel[3], double mag[3], double DCM[3][3]);
