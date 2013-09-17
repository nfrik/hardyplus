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
#include "plotter.h"
#include <unistd.h>


using namespace Eigen;
int main(int argc, const char * argv[])
{
    hardycpp *d = new hardycpp::hardycpp("/Users/nfrik/Documents/LAMMPS/work/couette_flow_git/flow.hpcvel.17.txt");
    
    unsigned long ts=time(NULL);
    d->run(1000, 10, 10, 1);
    cout<<"total time in sec: "<<time(NULL)-ts<<endl;
    
//      MatrixXd insiders;
//      MatrixXd outsiders;
//      d->getBodyHeadTail2Matrix(data, 1000, 0.1);
//      d->getInsideAtoms(data, insiders, true, 0.1, 0.3, 0.1, 0.3, -INFINITY, INFINITY);
//      d->getOutsideAtoms(data, outsiders, true, 0.1, 0.3, 0.1, 0.3, -INFINITY, INFINITY, 0.1, 0.1, 0.1);
//////
//////
//////    
//      d->plot(outsiders.col(1).data(), outsiders.col(2).data(), (int)outsiders.col(2).size());
//    Vector3d d(10,30,230);
//    Vector3d c(5,1,5);
////    c<<5,1,5;
//    MatrixXd m(3,3);
//    m<<1,2,3,
//       4,5,6,
//        7,8,9;
//    m=(m+m.transpose()).eval();
//    cout<<m<<endl;
// //   cout<<(((c.transpose()-m.row(0)).array().abs())-(5*c.transpose()-m.row(0)).array().abs()).sum()<<endl;
// //   sort(c.data(), c.data()+c.size());
////    d=c.array()-5;
//    c=VectorXd::Map(m.transpose().row(0).data(), 3);
//    cout<<c<<endl;
//

    delete d;
    
    return 0;
}
