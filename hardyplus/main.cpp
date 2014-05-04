//
//  main.cpp
//  hardyplus
//
//  Created by Nikolay Frik on 8/21/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#include <iostream>
#include <string>
#include <algorithm>
#include <math.h>
#include "hardycpp.h"
#include <Eigen/Dense>
//#include "plotter.h"
#include <unistd.h>


using namespace Eigen;

extern int compute(int argc, char * argv[]);
extern void test();


int main(int argc, char * argv[])
{
    
    //---------------------------------------------------------------------------------------------------------------
    unsigned long ts=time(NULL);
    //---------------------------------------------------------------------------------------------------------------
    
    char cwf[1024]; // current working file
    int nx,ny,nz,t0,t1;
    hardycpp *d;
    
    //build for run
    if (true) {
        
        if(getcwd(cwf, sizeof(cwf)) == NULL){
            perror("getcwd(...) error");
        }
        
        //
        //hardyplus 10 10 1 0 1000 /Users/nfrik/Documents/..../file
        
        if(argc<7){
            printf("%s %d","Error: not enough arguments \n ",argc);
            return 0;
        }
        
        sscanf(argv[1],"%d", &nx);
        sscanf(argv[2],"%d", &ny);
        sscanf(argv[3],"%d", &nz);
        sscanf(argv[4],"%d", &t0);
        sscanf(argv[5],"%d", &t1);
        
        if(argv[6][0]=='.')
            strcat(cwf, &argv[6][1]);
        else if(argv[6][0]=='/') //we assume user supplies full path
        {
            strcpy(cwf, &argv[6][0]);
        }
        else //we assume user doesn't supply full path
        {
            printf("Adding // char \n ");
            strcat(cwf,"/");
            strcat(cwf, &argv[6][0]);
        }
        
        
        //printf("Command received with arguments: %s %d %d %d %d %d %s \n",argv[0],nx,ny,nz,t0,t1,cwf);
        printf("Reading file: %s \n",cwf);
        
        d = new typename hardycpp::hardycpp(cwf);
        for (int i=t0; i<=t1; i+=1000) {
            cout<<"Run: "<<i<<endl;
            d->run(i, nx, ny, nz);
        }
        delete d;
        d=NULL;
    }//if<build for run>
    else{
        t0=0;
        t1=10000;
        nx=20;
        ny=10;
        nz=1;
        d = new typename hardycpp::hardycpp("/Users/nfrik/Documents/LAMMPS/work/01262014/ber2nh_02132014/procedural_test/flow.berendsen-nosehoover-wide.xtx");
        d->run(t0, nx, ny, nz);
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
    if (d!=NULL) {
        delete d;
    }
    
    //---------------------------------------------------------------------------------------------------------------
    cout<<"Task completed in: "<<time(NULL)-ts<<" seconds"<<endl;
    //---------------------------------------------------------------------------------------------------------------
    return 0;

    
    //test();
    //return compute(argc, argv);
}


//void test(){
//    
//    interactions *phys = new interactions();
//    Vector3d A(0,0,0);
//    Vector3d B(0,0,0);
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<0,1.0,0;
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<1.0,0,0;
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<-1.0,0.0,0;
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<0,-1.0,0;
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<0,2.5,0;
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<2.5,0,0;
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<0,-2.5,0;
//    cout << phys->LJForce(A, B, 2.2);
//    A<<0,0,0;
//    B<<-2.5,0,0;
//    cout << phys->LJForce(A, B, 2.2);
//}
//
//int compute(int argc, char * argv[]){
//    //---------------------------------------------------------------------------------------------------------------
//    unsigned long ts=time(NULL);
//    //---------------------------------------------------------------------------------------------------------------
//    
//    char cwf[1024]; // current working file
//    int nx,ny,nz,t0,t1;
//    hardycpp *d;
//    
//    //build for run
//    if (true) {
//        
//        if(getcwd(cwf, sizeof(cwf)) == NULL){
//            perror("getcwd(...) error");
//        }
//        
//        //
//        //hardyplus 10 10 1 0 1000 /Users/nfrik/Documents/..../file
//        
//        if(argc<7){
//            printf("%s %d","Error: not enough arguments \n ",argc);
//            return 0;
//        }
//        
//        sscanf(argv[1],"%d", &nx);
//        sscanf(argv[2],"%d", &ny);
//        sscanf(argv[3],"%d", &nz);
//        sscanf(argv[4],"%d", &t0);
//        sscanf(argv[5],"%d", &t1);
//        
//        if(argv[6][0]=='.')
//            strcat(cwf, &argv[6][1]);
//        else if(argv[6][0]=='/') //we assume user supplies full path
//        {
//            strcpy(cwf, &argv[6][0]);
//        }
//        else //we assume user doesn't supply full path
//        {
//            printf("Adding // char \n ");
//            strcat(cwf,"/");
//            strcat(cwf, &argv[6][0]);
//        }
//        
//        
//        //printf("Command received with arguments: %s %d %d %d %d %d %s \n",argv[0],nx,ny,nz,t0,t1,cwf);
//        printf("Reading file: %s \n",cwf);
//        
//        d = new hardycpp::hardycpp(cwf);
//        for (int i=t0; i<=t1; i+=1000) {
//            cout<<"Run: "<<i<<endl;
//            d->run(i, nx, ny, nz);
//        }
//        delete d;
//        d=NULL;
//    }//if<build for run>
//    else{
//        t0=0;
//        t1=10000;
//        nx=20;
//        ny=10;
//        nz=1;
//        d = new hardycpp::hardycpp("/Users/nfrik/Documents/LAMMPS/work/01262014/ber2nh_02132014/procedural_test/flow.berendsen-nosehoover-wide.xtx");
//        d->run(t0, nx, ny, nz);
//    }
//    
//    
//    
//    //      MatrixXd outsiders;
//    //      d->getBodyHeadTail2Matrix(data, 1000, 0.1);
//    //      d->getInsideAtoms(data, insiders, true, 0.1, 0.3, 0.1, 0.3, -INFINITY, INFINITY);
//    //      d->getOutsideAtoms(data, outsiders, true, 0.1, 0.3, 0.1, 0.3, -INFINITY, INFINITY, 0.1, 0.1, 0.1);
//    //////
//    //////
//    //////
//    //      d->plot(outsiders.col(1).data(), outsiders.col(2).data(), (int)outsiders.col(2).size());
//    //    Vector3d d(10,30,230);
//    //    Vector3d c(5,1,5);
//    ////    c<<5,1,5;
//    //    MatrixXd m(3,3);
//    //    m<<1,2,3,
//    //       4,5,6,
//    //        7,8,9;
//    //    double d=m(0,1)+m(0,1);
//    //    cout<<d<<endl;
//    //    m=(m+m.transpose()).eval();
//    //    cout<<m<<endl;
//    // //   cout<<(((c.transpose()-m.row(0)).array().abs())-(5*c.transpose()-m.row(0)).array().abs()).sum()<<endl;
//    // //   sort(c.data(), c.data()+c.size());
//    ////    d=c.array()-5;
//    //    c=VectorXd::Map(m.transpose().row(0).data(), 3);
//    //    cout<<c<<endl;
//    //
//    if (d!=NULL) {
//        delete d;
//    }
//
//    //---------------------------------------------------------------------------------------------------------------
//    cout<<"Task completed in: "<<time(NULL)-ts<<" seconds"<<endl;
//    //---------------------------------------------------------------------------------------------------------------        
//    return 0;
//}
