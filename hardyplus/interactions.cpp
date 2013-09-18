//
//  interactions.cpp
//  hardyplus
//
//  Created by Nikolay Frik on 9/10/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#include "interactions.h"
#include <math.h>
#include <float.h>

#define  DBLDBL_MIN 2*DBL_MIN


interactions::interactions(){
    
}

interactions::~interactions(){
    
}


// with default sigma=epsilon=1.0
Eigen::Vector3d interactions::LJForce(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc){
    return LJForce(A, B, rc, 1.0, 1.0);
}

// with default sigma=epsilon=1.0
double interactions::LJPotential(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc){
    return LJPotential(A, B, rc, 1.0, 1.0);
}

Eigen::Vector3d interactions::LJForce(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc, double sigma, double epsilon){

//    Eigen::Vector3d r=Eigen::sqrt(sum(((A-B).^2)'));
    double r=sqrt((A-B).array().square().sum());
    Eigen::Vector3d rr=(A-B)/r;
    Eigen::Vector3d f(0,0,0);
    if ((r<rc)&&(r>0))
      f=24*epsilon*rr*(2*pow(sigma/r, 12)/r-pow(sigma/r, 6)/r);
    return f;
}

double interactions::LJPotential(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc, double sigma, double epsilon){
    
    double r=sqrt((A-B).array().square().sum());
    double u=0;
    if ((r<rc)&&(r>0))
        u=4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
    return u;
}



// NEED TO BE TESTED!!!!!
//
// A is alway an outside particle
//
//
double interactions::getLambdaBond(const Eigen::Vector3d A, Eigen::Vector3d B, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax){

    if (abs(B(1)-A(1))<=DBLDBL_MIN) {
        //same y different x
        if (A(0)<xmin)
            return abs(B(0)-xmin)/abs(A(0)-B(0));
        else if (A(0)>xmax)
            return abs(B(0)-xmax)/abs(A(0)-B(0));
    }
    
    if (abs(B(0)-A(0))<=DBLDBL_MIN) {
        //same x different y
        if (A(1)<ymin)
            return abs(B(1)-ymin)/abs(A(1)-B(1));
        else if (A(1)>ymax)
            return abs(B(1)-ymax)/abs(A(1)-B(1));
    }
    
    double k=(B(1)-A(1))/(B(0)-A(0));
    double b=B(1)-k*B(0);
    
    double ym1=k*xmin+b;
    double xm1=xmin;
    double ym2=ymin;
    double xm2=(ymin-b)/k;
    double ym3=ymax;
    double xm3=(ymax-b)/k;
    double ym4=k*xmax+b;
    double xm4=xmax;
    int i,n;
    
    Eigen::MatrixXd M(4,3);
    M<<xm1,ym1,0,
       xm2,ym2,0,
       xm3,ym3,0,
       xm4,ym4,0;
    
    Eigen::Vector3d Mm(4);

    Mm(0)=((B.transpose()-M.row(0)).array().abs()+(A.transpose()-M.row(0)).array().abs()).sum();
    Mm(1)=((B.transpose()-M.row(1)).array().abs()+(A.transpose()-M.row(1)).array().abs()).sum();
    Mm(2)=((B.transpose()-M.row(2)).array().abs()+(A.transpose()-M.row(2)).array().abs()).sum();
    Mm(3)=((B.transpose()-M.row(3)).array().abs()+(A.transpose()-M.row(3)).array().abs()).sum();
    
    //Find minimal Mm
    for (n=i=0; i<4; i++)
        n=(Mm(i)<Mm(n))?i:n;
    
    return sqrt((B.transpose()-M.row(0)).array().square().sum())/sqrt((B-A).array().square().sum());
}
