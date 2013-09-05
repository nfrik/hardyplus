//
//  hardycpp.h
//  hardyplus
//
//  Created by Nikolay Frik on 8/21/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#ifndef __hardyplus__hardycpp__
#define __hardyplus__hardycpp__

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#endif /* defined(__hardyplus__hardycpp__) */

using namespace std;

struct simdata{
    int natoms; //number of atoms
    int timestep=1000; //timestep of output data
    double xmin,xmax,ymin,ymax,zmin,zmax;
};

class hardycpp{
public:
    //specify path and starting/ending t0/t1 parameters
    hardycpp(const char *path);
    ~ hardycpp();
    void test();
    void run(int time,int dargx, int dagy, int dargz);
    vector<int> findindxs(bool scaled, int time, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    void getBodyHeadTail2Matrix(Eigen::MatrixXd &m,int time, double rcx);
    void plot(const double *xData, const double *yData,int dataSize);
    int k;
private:
    simdata sdata;
    void readRawDataFromFile(const char *str);
    vector<vector<double>> wdata; //whole data
    vector<int> inatomsindx;
    vector<int> outatomsindx;
//    vector<vector<vector<double>>> inatoms;
//    vector<vector<vector<double>>> outatoms;
    
    int b;
};
