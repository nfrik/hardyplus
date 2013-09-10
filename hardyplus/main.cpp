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
    hardycpp *d = new hardycpp::hardycpp("/Users/nfrik/Documents/LAMMPS/work/couette_flow_git/flow.hpcvel.17.txt");
    
    MatrixXd data;
    MatrixXd insiders;
    MatrixXd outsiders;
    d->getBodyHeadTail2Matrix(data, 1000, 0.1);
//    d->getInsideAtoms(data, insiders, true, 0.1, 0.3, 0.1, 0.3, -INFINITY, INFINITY);
    d->getOutsideAtoms(data, outsiders, true, 0.1, 0.3, 0.1, 0.3, -INFINITY, INFINITY, 0.1, 0.1, 0.1);
//    
//
//    
    d->plot(outsiders.col(1).data(), outsiders.col(2).data(), (int)outsiders.col(2).size());

    
//    MatrixXd m(3,3);
//    m<<1,2,3,
//       4,5,6,
//        7,8,9;
//    
//    cout<<m.cols()<<endl;
//    
//    VectorXd v(3),b;
//    v << 1,2,3;
//
//    vector<vector<double>> s;
//    s.resize(2);
//    s[0].resize(2);
//    s[1].resize(2);
//    s[0][0]=1;
//    s[1][0]=2;
//    s[0][1]=3;
//    s[1][1]=4;
//    cout<<s.size()*s[0].size()<<endl;
//    Map<MatrixXd> mymat(&s[0][0],2,2);
//    cout<<mymat<<endl;
//    b=v;
//    b.transposeInPlace();
    //cout<<m.cwiseQuotient(m.transpose().eval())<<endl;
    
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
    
    
    //int n=5;
    
    //std::cout<<"Upper bound = "<<*it<<endl;
    

    
    
    delete d;
    
    return 0;
}
