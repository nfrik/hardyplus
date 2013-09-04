//
//  main.cpp
//  hardyplus
//
//  Created by Nikolay Frik on 8/21/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <math.h>
#include "hardycpp.h"
#include <Eigen/Dense>

using namespace Eigen;
int main(int argc, const char * argv[])
{
    //hardycpp *d = new hardycpp::hardycpp("/Users/nfrik/Documents/LAMMPS/work/couette_flow_git/flow.hpcvel.17.txt");
    
    MatrixXd m=MatrixXd::Constant(3, 3, 5);
    VectorXd v(3);
    v << 1,2,3;
    cout<<m*v<<endl;
//    vector<vector<int>> d;
//    
//    d.resize(10);
//    for (int i=0; i<10; i++) {
//        d[i].resize(3);
//        d[i][0]=i;
//        d[i][1]=i*10;
//        d[i][2]=i*100;
//    }
//
//    for (int i=0; i<10; i++) {
//        std::cout<<"d["<<i<<"][0]="<<d[i][0]<<"  d["<<i<<"][1]="<<d[i][1]<<"  d["<<i<<"][2]="<<d[i][2]<<endl;
//    }
    
    
//    std::cout<<"Lower bound"<<*low<<endl;
//    std::cout<<"Upper bound"<<*up<<endl;
    
    std::vector<vector<int>>::iterator row;
    
    //int n=5;
    
    //std::cout<<"Upper bound = "<<*it<<endl;
    

    
    
    //delete d;
    
    return 0;
}
