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
//---------------------------------------------------------------------------------------------------------------        
   unsigned long ts=time(NULL);
//---------------------------------------------------------------------------------------------------------------

    
   hardycpp *d = new hardycpp::hardycpp("/Users/nfrik/Documents/LAMMPS/work/couette_flow_git/flow.hpcvel.17.txt");
    
//    d->run(000, 10, 10, 1);
    
    for (int i=0; i<=200000; i+=1000) {
        cout<<"Frame: "<<i<<endl;
        d->run(i, 10, 10, 1);
    }
    
    
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
//    double d=m(0,1)+m(0,1);
//    cout<<d<<endl;
//    m=(m+m.transpose()).eval();
//    cout<<m<<endl;
// //   cout<<(((c.transpose()-m.row(0)).array().abs())-(5*c.transpose()-m.row(0)).array().abs()).sum()<<endl;
// //   sort(c.data(), c.data()+c.size());
////    d=c.array()-5;
//    c=VectorXd::Map(m.transpose().row(0).data(), 3);
//    cout<<c<<endl;
//

    delete d;
    
//---------------------------------------------------------------------------------------------------------------
    cout<<"Task completed in: "<<time(NULL)-ts<<" seconds"<<endl;
//---------------------------------------------------------------------------------------------------------------        
    return 0;
}
