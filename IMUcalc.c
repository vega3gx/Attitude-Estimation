#include <math.h>
#include "BOARD.h"
#include "IMUcalc.h"
#include <stdio.h>
#include "MatrixMath.h"

#define KPA 10
#define KIA 5
#define KPM 10
#define KIM 5
#define GRAV_INERTIAL {0,0,1}
#define MAG_INERTIAL {0.4772,0.1128,0.8715}

double psiFromDCM(double DCMMatrix[3][3]){
    double Euler[3];
    EulFromDCM(DCMMatrix,Euler);
    return(180.0*Euler[0]/M_PI);
}

double phiFromDCM(double DCMMatrix[3][3]){
    double Euler[3];
    EulFromDCM(DCMMatrix,Euler);
    return(180.0*Euler[2]/M_PI);
}

double thetaFromDCM(double DCMMatrix[3][3]){
    double Euler[3];
    EulFromDCM(DCMMatrix,Euler);
    return(180.0*Euler[1]/M_PI);
}

void EulFromDCM(double DCMMatrix[3][3],double Euler[3]){
    if(DCMMatrix[0][2]!=1 && DCMMatrix[0][2]!=-1){
        Euler[1]=-1.0*asin(DCMMatrix[0][2]);
        Euler[0]=atan2(DCMMatrix[0][1],DCMMatrix[1][1]);
        Euler[2]=atan2(DCMMatrix[1][2],DCMMatrix[2][2]);      
    }else{
        Euler[2]=0;
        if(DCMMatrix[0][2]==-1){
            Euler[1]=M_PI/2.0;
            Euler[0]=atan2(DCMMatrix[2][1],DCMMatrix[2][0]);
        }else{
            Euler[1]=-1.0*M_PI/2.0;
            Euler[0]=atan2(-1.0*DCMMatrix[2][1],DCMMatrix[2][0]);
        }
    }
}

double sinc(double x){
    if(x<0.5 && x>-0.5){return(1-(M_PI*M_PI*x*x/6)+(M_PI*M_PI*M_PI*M_PI*x*x*x*x/120));}
    else{return(sin(M_PI*x)/(M_PI*x));}
}

void Rexp(double w[3], double ans[3][3]){
    double wnorm=norm(w);
    double rx[3][3] = {
        {0,-1.0*w[2],w[1]},
        {w[2],0,-1.0*w[0]},
        {-1.0*w[1],w[0],0}
    };
    double I[3][3] ={
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };
    double s=sinc(wnorm/2/M_PI);
    double c=cos(wnorm/2);
    double rx2[3][3];
    double rx3[3][3];
    double rx4[3][3];
    double rx5[3][3];
    MatrixMultiply(rx,rx,rx2);
    MatrixScalarMultiply(s*s*0.5, rx2, rx3);
    MatrixScalarMultiply(s*c,rx,rx4);
    MatrixAdd(rx3,rx4,rx5);
    MatrixAdd(rx5,I,ans);
}

double norm(double x[3]){
    return(sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)));
}

void IntegrateOpenLoopMatExp(double Rminus[3][3], double gyro[3], double deltaT, double Rplus[3][3]){
    double rexp[3][3];
    double scaledGyro[3];
    int i;
    for(i=0; i<3; i++){scaledGyro[i]=gyro[i]*deltaT;}
    Rexp(scaledGyro,rexp);
    MatrixMultiply(Rminus,rexp,Rplus);
}

void IntegrateOpenLoopForwardInt(double Rminus[3][3], double gyro[3], double deltaT, double Rplus[3][3]){
    double scaledGyro[3];
    int i;
    for(i=0; i<3; i++){scaledGyro[i]=gyro[i]*deltaT;}
    double rx[3][3] = {
        {0,-1.0*scaledGyro[2],scaledGyro[1]},
        {scaledGyro[2],0,-1.0*scaledGyro[0]},
        {-1.0*scaledGyro[1],scaledGyro[0],0}
    };
    double rx2[3][3];
    MatrixMultiply(Rminus,rx,rx2);
    MatrixAdd(Rminus,rx2,Rplus);
}

void IntegrateClosedLoop(double Rminus[3][3], double Bminus[3], double gyros[3], double mags[3], double accel[3], double dT, double Rplus[3][3], double Bplus[3]){
    double aNorm=norm(accel);
    double magNorm=norm(mags);
    int i;
    for(i=0;i<3;i++){
        accel[i]=accel[i]/aNorm;
        mags[i]=mags[i]/magNorm;
    }
    double gyroInWithBias[3];
    double negBias[3];
    VectorScale(Bminus,-1.0,negBias);
    VectorAdd(gyros,negBias,gyroInWithBias);
    double rxa[3][3] = {
        {0,-1.0*accel[2],accel[1]},
        {accel[2],0,-1.0*accel[0]},
        {-1.0*accel[1],accel[0],0}
    };
    double rxm[3][3] = {
        {0,-1.0*mags[2],mags[1]},
        {mags[2],0,-1.0*mags[0]},
        {-1.0*mags[1],mags[0],0}
    };
    double wmeas_a[3];
    double wmeas_m[3];
    double aInert[3]=GRAV_INERTIAL;
    double mInert[3]=MAG_INERTIAL;
    double aBody[3];
    double mBody[3];
    double invRmin[3][3];
    MatrixTranspose(Rminus,invRmin);
    MatrixVectorProduct(invRmin,aInert,aBody);
    MatrixVectorProduct(invRmin,mInert,mBody);
    MatrixVectorProduct(rxa,aBody,wmeas_a);
    MatrixVectorProduct(rxm,mBody,wmeas_m);
    double gyroWFeedback[3];
    double feedback[3];
    double Afeedback[3];
    double Mfeedback[2];
    VectorScale(wmeas_a,KPA,Afeedback);
    VectorScale(wmeas_m,KPM,Mfeedback);
    VectorAdd(Afeedback,Mfeedback,feedback);
    VectorAdd(feedback,gyroInWithBias,gyroWFeedback);
    double bdot[3];
    double Aifeedback[3];
    double Mifeedback[3];
    VectorScale(wmeas_a,-1*KIA,Aifeedback);
    VectorScale(wmeas_m,-1*KIM,Mifeedback);
    VectorAdd(Aifeedback,Mifeedback,bdot);
    double rexp[3][3];
    double scaleGyro[3];
    VectorScale(gyroWFeedback,dT,scaleGyro);
    Rexp(scaleGyro,rexp);
    MatrixMultiply(Rminus,rexp, Rplus);
    double sbdot[3];
    VectorScale(bdot,dT,sbdot);
    VectorAdd(Bminus,sbdot,Bplus);
}

void DCMfromTriad(double accel[3], double mag[3], double DCM[3][3]){
    double aNorm=norm(accel);
    double magNorm=norm(mag);
    int i;
    for(i=0;i<3;i++){
        accel[i]=accel[i]/aNorm;
        mag[i]=mag[i]/magNorm;
    }
    double rxa[3][3] = {
        {0,-1.0*accel[2],accel[1]},
        {accel[2],0,-1.0*accel[0]},
        {-1.0*accel[1],accel[0],0}
    };
    double rxm[3][3] = {
        {0,-1.0*mag[2],mag[1]},
        {mag[2],0,-1.0*mag[0]},
        {-1.0*mag[1],mag[0],0}
    };
    double m[3];
    double M[3];
    double sm[3];
    double sM[3];
    MatrixVectorProduct(rxm,accel,m);
    MatrixVectorProduct(rxa,mag,M);
    VectorScale(m,1/norm(m),sm);
    VectorScale(M,1/norm(M),sM);
    double aInert[3]=GRAV_INERTIAL;
    double mInert[3]=MAG_INERTIAL;
    double rxai[3][3] = {
        {0,-1.0*aInert[2],aInert[1]},
        {aInert[2],0,-1.0*aInert[0]},
        {-1.0*aInert[1],aInert[0],0}
    };
    double rxmi[3][3] = {
        {0,-1.0*mInert[2],mInert[1]},
        {mInert[2],0,-1.0*mInert[0]},
        {-1.0*mInert[1],mInert[0],0}
    };
    double mxi[3];
    double Mxi[3];
    MatrixVectorProduct(rxmi,M,Mxi);
    MatrixVectorProduct(rxai,m,mxi);
    double mul1[3][3];
    double mul2[3][3];
    for(i=0;i<3;i++){
        mul1[i][0]=mInert[i];
        mul2[i][0]=mag[i];
        mul1[i][1]=M[i];
        mul2[i][1]=m[i];
        mul1[i][2]=Mxi[i];
        mul2[i][2]=mxi[i];
    }
    MatrixMultiply(mul1,mul2,DCM);
}
