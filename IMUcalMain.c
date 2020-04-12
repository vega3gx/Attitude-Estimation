/*
 * File:   IMUcalMain.c
 * Author: aamwood
 *
 * Created on February 18, 2020, 2:43 PM
 */


#include "Ascii.h"
#include "BOARD.h"
#include "I2C.h"
#include "MPU9250.h"
#include "serial.h"
#include "timers.h"
#include "IMUcalc.h"
#include "MatrixMath.h"
#include <stdio.h>

#define SAMPLE_PERIOD 50
#define BIAS_SAMPLES 20000
#define GYRO2RAD_PER_SEC 0.000266462
#define ACCEL_SCALE 16834
#define MAG_SCALE 1666.66
//#define EULTEST
//#define SINCTEST
//#define REXPTEST
//#define OLGYRO
//#define DATAOUT
#define CLGYRO
//#define CLTEST
//#define MATMATH
//#define TRIAD

double magMisAlign[3][3] ={
    {0.9999,-0.0007,0.0128},
    {0.0006,1.0000,0.0068},
    {-0.0128,-0.0068,0.9999}
};

double AccelA[3][3] ={
    {1.0059,.00007583,0.0133},
    {0.0049,1.0014,-0.0080},
    {0.0014,-0.0128,0.9912}
};

double AccelB[3] = {-0.0510,-0.0161,-0.0341};

double MagA[3][3] = {
    {0.1530,-0.0015,-0.0093},
    {-0.0014,0.1498,0.0026},
    {-0.0106,-0.0003414,0.1354}
};

double MagB[3] = {-0.8452,-0.3857,0.5006};

int main(void) {
    BOARD_Init();
    printf("IMU Data \n");
    MPU9250_Init();
    TIMERS_Init();
    OledInit();
#ifdef DATAOUT
    while(1){
        int ax,ay,az,gx,gy,gz,mx,my,mz;
            while(TIMERS_GetMilliSeconds()%SAMPLE_PERIOD);
            gx=MPU9250_ReadGyroX(); 
            gy=MPU9250_ReadGyroY();
            gz=MPU9250_ReadGyroZ();
            ax=MPU9250_ReadAccelX(); 
            ay=MPU9250_ReadAccelY();
            az=MPU9250_ReadAccelZ();
            mx=MPU9250_ReadMagX(); 
            my=MPU9250_ReadMagY();
            mz=MPU9250_ReadMagZ();
        printf("%d,%d,%d,%d,%d,%d,%d,%d,%d\r\n",ax,ay,az,gx,gy,gz,mx,my,mz);
        //printf("%d,%d,%d\r\n",mx,my,mz);
    }
#endif
#ifdef EULTEST
        double DCM[3][3] ={
                            {-0.2539,  0.0870, -0.9633},
                            { 0.9672, -0.6065, -0.2534},
                            { 0.6063, -0.9961, -0.0884}
        };
        printf("psi = \n%f \n",psiFromDCM(DCM));
//      should be 161.09
        printf("theta = %f \n\n",thetaFromDCM(DCM));
//      should be 74.4335
        printf("phi = %f \n\n",phiFromDCM(DCM));
//      should be -109.24
       while(1);
#endif
#ifdef SINCTEST
       printf("sinc(0.3) = %f \n\n",sinc(0.3));
       printf("sinc(0.7) = %f \n\n",sinc(0.7));
       printf("sinc(1.5) = %f \n\n",sinc(1.5));
       while(1);
#endif
#ifdef REXPTEST
       double w[3]={100,50,-200};
       double rexp[3][3];
       Rexp(w,rexp);
       printf("result = \r\n");
       MatrixPrint(rexp);
       while(1);
#endif
#ifdef OLGYRO   
       printf("olgyro\n\r");
    int i;
    int bias[3] = {0,0,0};
    double dt = 1.0/SAMPLE_PERIOD;
    for(i=0;i<BIAS_SAMPLES;i++){
        bias[1]+=MPU9250_ReadGyroX();
        bias[2]+=MPU9250_ReadGyroY();
        bias[3]+=MPU9250_ReadGyroZ();
    }
    double DCM[3][3] = {
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };
    double out[3][3]= {
        {1,0,0},
        {0,1,0},
        {0,0,1}};
    for(i=0;i<3;i++){bias[i]/=BIAS_SAMPLES;}
    while(1){
        while(TIMERS_GetMilliSeconds()%SAMPLE_PERIOD);
//        double DCM[3][3] = {
//        {1,0,0},
//        {0,1,0},
//        {0,0,1}};
        double gyro[3]={10,10,10};
        gyro[0]=MPU9250_ReadGyroX()-bias[0]; 
        gyro[1]=MPU9250_ReadGyroY()-bias[1];
        gyro[2]=MPU9250_ReadGyroZ()-bias[2];
        for(i=0;i<3;i++){gyro[i]*=GYRO2RAD_PER_SEC;}
        IntegrateOpenLoopMatExp(DCM,gyro,dt,out);
        MatrixCopy(out,DCM);
        double psi = psiFromDCM(out);
        double theta = thetaFromDCM(out);
        double phi = phiFromDCM(out);
        //printf("Mat dcm = \r\n");
        //MatrixPrint(out);
        char oled[50];
        sprintf(oled,"psi = %3.2f \ntheta = 3.2%f \nphi = %3.2f \n",psi,theta,phi);
        OledClear(0);
        OledDrawString(oled);
        OledUpdate();
       }
#endif
#ifdef CLGYRO
    printf("clgyro\r\n");
    int i;
    double dt = 1.0/SAMPLE_PERIOD;
    double DCM[3][3] = {
        {1,0,0},
        {0,1,0},
        {0,0,1}
    };
    double out[3][3]= {
        {1,0,0},
        {0,1,0},
        {0,0,1}};
    double bias[3]={0,0,0};
    double bout[3]={0,0,0};
    while(1){
        while(TIMERS_GetMilliSeconds()%SAMPLE_PERIOD);
//        double DCM[3][3] = {
//        {1,0,0},
//        {0,1,0},
//        {0,0,1}};
        double gyro[3]={10,10,10};
        double accel[3];
        double mag[3];
        gyro[0]=MPU9250_ReadGyroX(); 
        gyro[1]=MPU9250_ReadGyroY();
        gyro[2]=MPU9250_ReadGyroZ();
        accel[0]=MPU9250_ReadAccelX(); 
        accel[1]=MPU9250_ReadAccelY();
        accel[2]=MPU9250_ReadAccelZ();
        mag[0]=MPU9250_ReadMagX(); 
        mag[1]=MPU9250_ReadMagY();
        mag[2]=MPU9250_ReadMagZ();
        
        
        for(i=0;i<3;i++){
            gyro[i]*=GYRO2RAD_PER_SEC;
            accel[i]/=ACCEL_SCALE;
            mag[i]*=MAG_SCALE;
        }
        
        double sacc[3];
        double smag[3];
        MatrixVectorProduct(AccelA,accel,sacc);
        MatrixVectorProduct(MagA,mag,smag);
        
        double sbacc[3];
        double sbmag[3];
        VectorAdd(sacc,AccelB,sbacc);
        VectorAdd(smag,MagB,sbmag);
        
        double magAligned[3];
//        printf("\r\n\r\n");
//        printf("bias = ");
//        VectorPrint(bias);
//        printf("gyro = ");
//        VectorPrint(gyro);
        printf("sbmag = ");
        VectorPrint(sbmag);
//        printf("sbacc = ");
//        VectorPrint(sbacc);
//        printf("DCMin = \r\n");
//        MatrixPrint(DCM);
        IntegrateClosedLoop(DCM,bias,gyro,sbmag,sbacc,dt,out,bout);
        MatrixCopy(out,DCM);
        VectorCopy(bout,bias);
        double psi = psiFromDCM(out);
        double theta = thetaFromDCM(out);
        double phi = phiFromDCM(out);
        //printf("Mat dcmOut = \r\n");
       // MatrixPrint(out);
       // printf("Bias Out = \r\n");
        //VectorPrint(bout);
        char oled[50];
        sprintf(oled,"psi = %3.2f \ntheta = %3.2f \nphi = %3.2f \n",psi,theta,phi);
        OledClear(0);
        OledDrawString(oled);
        OledUpdate();
       }
#endif
#ifdef CLTEST
    double dt = .02;
    double bias[3] = {0,0,0};
    double bplus[3];
    double rplus[3][3];
    double DCMin[3][3] = {
        {.9990,.0113,-.0422},
        {-.0117,.9999,-.0094},
        {.0421,.0099,.9991}
    };
    double gyro[3] = {-.7346, .1037, .2635};
    double sbmag[3] = {-.8452, -.3857, .5006};
    double sbacc[3] = {.5645, .1687, .8585};
    IntegrateClosedLoop(DCMin,bias,gyro,sbmag,sbacc,dt,rplus,bplus);
    printf("Mat dcmOut = \r\n");
    MatrixPrint(rplus);
    printf("Bias Out = \r\n");
    VectorPrint(bplus);
    while(1);
#endif
#ifdef MATMATH
    double mat[3][3] = {
        {.9990,.0113,-.0422},
        {-.0117,.9999,-.0094},
        {.0421,.0099,.9991}
    };
    double vect[3] = {0,0,1};
    double out[3];
    MatrixVectorProduct(mat,vect,out);
    VectorPrint(out);
#endif
#ifdef TRIAD
    printf("triad\n\r");
    int i;
    double dt = 1.0/SAMPLE_PERIOD;
    double out[3][3]= {
        {1,0,0},
        {0,1,0},
        {0,0,1}};
    while(1){
        while(TIMERS_GetMilliSeconds()%SAMPLE_PERIOD);
//        double DCM[3][3] = {
//        {1,0,0},
//        {0,1,0},
//        {0,0,1}};
        double mag[3];
        double accel[3];
        accel[0]=MPU9250_ReadAccelX(); 
        accel[1]=MPU9250_ReadAccelY();
        accel[2]=MPU9250_ReadAccelZ();
        mag[0]=MPU9250_ReadMagX(); 
        mag[1]=MPU9250_ReadMagY();
        mag[2]=MPU9250_ReadMagZ();
        
        for(i=0;i<3;i++){
            accel[i]/=ACCEL_SCALE;
            mag[i]*=MAG_SCALE;
        }
        
        double sacc[3];
        double smag[3];
        MatrixVectorProduct(AccelA,accel,sacc);
        MatrixVectorProduct(MagA,mag,smag);
        
        double sbacc[3];
        double sbmag[3];
        VectorAdd(sacc,AccelB,sbacc);
        VectorAdd(smag,MagB,sbmag);
        
        DCMfromTriad(sbacc,sbmag,out);
        double psi = psiFromDCM(out);
        double theta = thetaFromDCM(out);
        double phi = phiFromDCM(out);
        //printf("Mat dcm = \r\n");
        //MatrixPrint(out);
        char oled[50];
        sprintf(oled,"psi = %3.2f \n theta = %3.2f \n phi = %3.2f \n",psi,theta,phi);
        OledClear(0);
        OledDrawString(oled);
        OledUpdate();
       }
#endif
    //}
}
