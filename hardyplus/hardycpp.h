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
    
    //move below functions to private after testing
    vector<int> findindxs(bool scaled, int time, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    void        findatoms(Eigen::MatrixXd &atoms,bool scaled, int time, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    
    void getBodyHeadTail2Matrix(Eigen::MatrixXd &m,int time, double rcx);

    //All address variables are output ones
    void getInsideAtoms(const Eigen::MatrixXd &data, Eigen::MatrixXd &atoms, bool scaled, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    void getOutsideAtoms(const Eigen::MatrixXd &data, Eigen::MatrixXd &atoms, bool scaled, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rcx, double rcy, double rcz);
    void neighborList(const Eigen::MatrixXd &InsidersIn, const Eigen::MatrixXd &OutsidersIn, Eigen::MatrixXd &phiOut, Eigen::MatrixXd &FxOut, Eigen::MatrixXd &FyOut, Eigen::MatrixXd &FzOut,
                      Eigen::MatrixXd &xijOut, Eigen::MatrixXd &yijOut, Eigen::MatrixXd &zijOut, Eigen::MatrixXd &lamOut);
    void stresskinetic(const Eigen::MatrixXd &InsidersIn, double avvxIn, double avvyIn, double avvzIn, double volIn, Eigen::MatrixXd &SkOut);
    void stresspotential(const Eigen::MatrixXd &FxIn, const Eigen::MatrixXd &FyIn, const Eigen::MatrixXd &FzIn,
                   const Eigen::MatrixXd &xijIn, const Eigen::MatrixXd &yijIn, const Eigen::MatrixXd &zijIn, const Eigen::MatrixXd &lamIn, double volIn, Eigen::MatrixXd &SpOut);
    
    void plot(const double *xData, const double *yData,int dataSize);
    void printMat2File(const Eigen::MatrixXd &m, string filename);
    
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
