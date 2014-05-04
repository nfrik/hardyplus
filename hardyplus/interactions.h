//
//  interactions.h
//  hardyplus
//
//  Created by Nikolay Frik on 9/10/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#ifndef __hardyplus__interactions__
#define __hardyplus__interactions__

#include <iostream>
#include <Eigen/Dense>
#include <vector>

#endif /* defined(__hardyplus__interactions__) */


using namespace std;
class interactions{
public:
    interactions();
    ~interactions();
    Eigen::Vector3d LJForce(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc);
    double LJPotential(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc);
    Eigen::Vector3d LJForce(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc, double sigma, double epsilon);
    double LJPotential(const Eigen::Vector3d A, const Eigen::Vector3d B, double rc, double sigma, double epsilon);
    double getLambdaBond(const Eigen::Vector3d A, Eigen::Vector3d B, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    bool isPMemberOfRectAB(double x, double y, double x0, double y0, double x1, double y1);
    double lam(double x0,double y0,double x1, double y1,double x2, double y2);
private:
    //add private members i.e. saving cutoff  variables
};
