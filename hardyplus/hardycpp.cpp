//
//  hardycpp.cpp
//  hardyplus
//
//  Created by Nikolay Frik on 8/21/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#include "hardycpp.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;

hardycpp::hardycpp(const char *path){
    
    
    this->readRawDataFromFile(path);//read file to wdata
    
    //this->run(10000,10,10,1);//run 10000'th frame with 10x10x1 mesh grid
    
}

hardycpp::~hardycpp(){
    
}

void hardycpp::run(int time, int dargx, int dargy, int dargz){
    cout<<"Natoms: "<<sdata.natoms<<"\nxmin xmax: "<<sdata.xmin<<" "<<sdata.xmax<<"\nymin ymax: "<<sdata.ymin<<" "<<sdata.ymax<<"\nzmin zmax: "<<sdata.zmin<<" "<<sdata.zmax<<endl;
    
    //wee need to compensate atom coordinates because of removal of walls
    double yclo=0.06675+0.00001;
    double ychi=0.9166-0.00001;
    double dely=ychi-yclo;
    
    
    //Recalculate rc (potential cutoff) from lj units and true box dimension
    double rc=1.12246;//potential cutoff distance
    double rcx=rc/(sdata.xmax-sdata.xmin);
    double rcy=rc/(sdata.ymax-sdata.ymin);
    double rcz=rc/(sdata.zmax-sdata.zmin);
    
    double sxlo,sylo,szlo,sxhi,syhi,szhi,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi,dsx,dsy,dsz,vol;
    
    //data preprocessing
    //pad data with tail and head in x direction
    //1***************2%     %2***************1%
    //1***************2%     %2***************1%
    //1***************2% --> %2***************1%
    //1***************2%     %2***************1%
    //1***************2%     %2***************1%
    //---------------------------
    
    
    MatrixXd data;
    getBodyHeadTail2Matrix(data, time, rcx);

    for (int i=1; i<=dargx; i++) {
        for (int j=1; j<=dargy; j++) {
            for (int k=1; k<=dargz; k++) {
                sxlo=(j-1)/dargx;
                sylo=(i-1)*dely/dargy+yclo;
                szlo=(k-1)/dargz;
                sxhi=(j)/dargx;
                syhi=(i)*dely/dargy+yclo;
                szhi=(k)/dargz;
                
                dx=(sdata.xmax-sdata.xmin)*(sxhi-sxlo);
                dy=(sdata.ymax-sdata.ymin)*(syhi-sylo);
                dz=(sdata.zmax-sdata.zmin)*(szhi-szlo);
                
//                dsx=(sxhi-sxlo);
//                dsy=(syhi-sylo);
//                dsz=(szhi-szlo);
                
                //true box dimensions for search
                xlo=sdata.xmin+dx*(j-1);
                xhi=sdata.xmin+dx*j;
                ylo=sdata.ymin+dy*(i-1)+(sdata.ymax-sdata.ymin)*yclo;
                yhi=sdata.ymin+dy*i+(sdata.ymax-sdata.ymin)*yclo;
                zlo=sdata.zmin+dz*(k-1);
                zhi=sdata.zmin+dz*k;
                
                vol=dx*dy*dz;
                
                vector<int> inatoms=findindxs(true, time, sxlo, sxhi, sylo, syhi, szlo, szhi);
                vector<int> outatoms=findindxs(true, time, sxlo-rcx, sxhi+rcx, sylo-rcy, syhi+rcy, szlo-rcz, szlo+rcz);
                
                
            }
        }
    }
    
//    double *xData=(double*)malloc(sizeof(double)*indxs.size());
//    double *yData=(double*)malloc(sizeof(double)*indxs.size());
//    for (int i=0; i<indxs.size(); i++) {
//        xData[i]=wdata[indxs[i]][7];
//        yData[i]=wdata[indxs[i]][8];
//        cout<<"Number i="<<i<<" -> "<<indxs[i]<<endl;
//    }
//    plot(xData, yData, (int)indxs.size());
//    delete []yData;
//    delete []xData;
}


void hardycpp::readRawDataFromFile(const char *str){
    
    ifstream readFileName;
    ///readFileName.exceptions(ifstream::badbit|ifstream::failbit|ifstream::eofbit);
    //char line[128];
    string line;
    
    //determine number of lines in file
    int flines=0;
    readFileName.open(str);
    if (readFileName.is_open())
        while (getline(readFileName, line))
            flines++;
    
    //std::count(std::istreambuf_iterator<char>(readFileName), std::istreambuf_iterator<char>(), '\n');

    //resize our matrix for
    this->wdata.resize(flines); // set height
    for(int i=0;i<flines;i++)
        this->wdata[i].resize(12); // set width
    
    readFileName.clear();
    readFileName.seekg(0,ios::beg);
    
    int i;
    if (readFileName.is_open()) {
        //read preamble
        for (i=0; i<7; i++){
            getline(readFileName, line);
            switch (i) {
                case 1:
                    sscanf(line.c_str(), "%d", &sdata.natoms);
                    break;
                case 3:
                    sscanf(line.c_str(), "%lf %lf", &sdata.xmin, &sdata.xmax);
                    break;
                case 4:
                    sscanf(line.c_str(), "%lf %lf", &sdata.ymin, &sdata.ymax);
                    break;
                case 5:
                    sscanf(line.c_str(), "%lf %lf", &sdata.zmin, &sdata.zmax);
                    break;
                default:
                    break;
            }
        }
        i=0;
        while (getline(readFileName, line)) {
              sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&wdata[i][0],&wdata[i][1],&wdata[i][2],&wdata[i][3],&wdata[i][4],&wdata[i][5],&wdata[i][6],&wdata[i][7],&wdata[i][8],&wdata[i][9],&wdata[i][10],&wdata[i][11]);
              i++;
        }
    }
    readFileName.close();
//    for (int i=10000;i<20000;i+=1000)
// //       cout<<wdata.size()<<endl;
//        cout<<wdata[i][0]<<" "<<wdata[i][1]<<" "<<wdata[i][2]<<" "<<wdata[i][3]<<" "<<wdata[i][4]<<" "<<wdata[i][5]<<" "<<wdata[i][6]<<" "<<wdata[i][7]<<" "<<wdata[i][8]<<" "<<wdata[i][9]<<" "<<wdata[i][10]<<" "<<wdata[i][11]<<endl;
}

// scaled - whether to search in scaled coordinates or real coordinates columns
// time - should be equal to real timestep you want to refer to

vector<int>  hardycpp::findindxs(bool scaled, int time,double xmin, double xmax, double ymin, double ymax, double zmin, double zmax){
    vector<int> indexes;
    int i; //timestep by default is scale of 1000
    const int start=time/sdata.timestep*sdata.natoms;
    const int end=(time/sdata.timestep+1)*sdata.natoms;
    
    
    if (scaled)//select what coordinates we need
    for (i=start; i<end; i++) {
        if ((xmin<=wdata[i][1])&&(wdata[i][1]<=xmax)&&(ymin<=wdata[i][2])&&(wdata[i][2]<=ymax)&&(zmin<=wdata[i][3])&&(wdata[i][3]<=zmax)) {
            indexes.push_back(i);
        }
    }
    else
    for (i=start; i<end; i++) {
        if ((xmin<=wdata[i][7])&&(wdata[i][7]<=xmax)&&(ymin<=wdata[i][8])&&(wdata[i][8]<=ymax)&&(zmin<=wdata[i][9])&&(wdata[i][9]<=zmax)) {
            indexes.push_back(i);            
        }
    }
    
    return indexes;
}

//Function assembles particles within cutoff range from left and right borders and returns it through matrix argument
//Matrix returned through input argument m,
//rcx - is a cutoff distance of lj potential
void hardycpp::getBodyHeadTail2Matrix(Eigen::MatrixXd &m,int time,double rcx){
    vector<int> indexHead;
    vector<int> indexTail;
    int i; //timestep by default is scale of 1000
    const int start=time/sdata.timestep*sdata.natoms;
    const int end=(time/sdata.timestep+1)*sdata.natoms;
    
//    vector<int> headIndxs=findindxs(true, time, 1.0-rcx, 1.0, -INFINITY, INFINITY, -INFINITY, INFINITY);
//    vector<int> tailIndxs=findindxs(true, time, 0.0, rcx, -INFINITY, INFINITY, -INFINITY, INFINITY);
    
    for (i=start; i<end; i++) {
        if ((1.0-rcx<=wdata[i][1])&&(wdata[i][1]<=1.0))
            indexHead.push_back(i);
        else if((0.0<=wdata[i][1])&&(wdata[i][1]<=rcx))
            indexTail.push_back(i);
    }
    
    unsigned long size=(indexHead.size()+indexTail.size()+end-start);
    
    m.resize(size, 11);
    
    int j=0;
    
//    head(:,2)=head(:,2)-1;
//    head(:,8)=head(:,8)-(sdata(1,2)-sdata(1,1));
//    tail(:,2)=tail(:,2)+1;
//    tail(:,8)=tail(:,8)+(sdata(1,2)-sdata(1,1));
    
    VectorXd headtailFix(11);
    headtailFix<<0,1,0,0,0,0,0,(sdata.xmax-sdata.xmin),0,0,0;
    
    //piecewize initialization
    for (unsigned long i=start; i<end; i++) {
        m.row(j)=VectorXd::Map(&wdata[i][0], 11);
        j++;
    }
    
    for(unsigned long i=0;i<indexHead.size();i++){
        m.row(j)=VectorXd::Map(&wdata[indexHead[i]][0], 11)-headtailFix;
        j++;
    }
    
    for(unsigned long i=0;i<indexTail.size();i++){
        m.row(j)=VectorXd::Map(&wdata[indexTail[i]][0], 11)+headtailFix;
        j++;
    }
    
}

void hardycpp::getInsideAtoms(const Eigen::MatrixXd &data, Eigen::MatrixXd &atoms, bool scaled, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax){

    vector <int> indexes;
    long end=data.rows();
    int i; //timestep by default is scale of 1000

    data(0,1);
    
    if (scaled)//select what coordinates we need
        for (i=0; i<end; i++) {
            if ((xmin<=data(i,1))&&(data(i,1)<=xmax)&&(ymin<=data(i,2))&&(data(i,2)<=ymax)&&(zmin<=data(i,3))&&(data(i,3)<=zmax)) {
                indexes.push_back(i);
            }
        }
    else
        for (i=0; i<end; i++) {
            if ((xmin<=data(i,7))&&(data(i,7)<=xmax)&&(ymin<=data(i,8))&&(data(i,8)<=ymax)&&(zmin<=data(i,9))&&(data(i,9)<=zmax)) {
                indexes.push_back(i);
            }
        }
    
    end=indexes.size();
    atoms.resize(end, 11);
    for (i=0; i<end; i++) {
//        for (int j=0; j<11; j++) {
//            atoms(i,j)=data(indexes[i],j);
//        }
        atoms.row(i)=data.row(indexes[i]);
        cout<<atoms.row(i)<<endl;
    }
    
}

void hardycpp::getOutsideAtoms(const Eigen::MatrixXd &data, Eigen::MatrixXd &atoms, bool scaled, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rcx, double rcy, double rcz){
    
    vector <int> indexes;
    long end=data.rows();
    int i; //timestep by default is scale of 1000
    
    //data(0,1);
    
    if (scaled)//select what coordinates we need
        for (i=0; i<end; i++) {
            ///Not a 3D case!!!!
            if (((xmin-rcx<=data(i,1))&&(data(i,1)<=xmin)&&(ymin-rcy<=data(i,2))&&(data(i,2)<=ymax+rcy))||
                ((xmax<=data(i,1))&&(data(i,1)<=xmax+rcx)&&(ymin-rcy<=data(i,2))&&(data(i,2)<=ymax+rcy))||
                ((xmin<data(i,1))&&(data(i,1)<xmax)&&(ymin-rcy<=data(i,2))&&(data(i,2)<=ymin))||
                ((xmin<data(i,1))&&(data(i,1)<xmax)&&(ymax<=data(i,2))&&(data(i,2)<=ymax+rcy))) {
                
                indexes.push_back(i);
            }
            
            
        }
    else{
        cout<<"Error in getOutsideAtoms(...) can't locate atoms in unscaled region"<<endl;
        exit(1);
    }
    
    end=indexes.size();
    atoms.resize(end, 11);
    for (i=0; i<end; i++) {
//        for (int j=0; j<11; j++) {
//            atoms(i,j)=data(indexes[i],j);
//        }
        atoms.row(i)=data.row(indexes[i]);
        cout<<atoms.row(i)<<endl;
    }
    
}

void hardycpp::neighborList(const Eigen::MatrixXd &InsidersIn, const Eigen::MatrixXd &OutsidersIn, Eigen::MatrixXd &phiOut, Eigen::MatrixXd &FxOut, Eigen::MatrixXd &FyOut, Eigen::MatrixXd &FzOut,
                  Eigen::MatrixXd &xijOut, Eigen::MatrixXd &yijOut, Eigen::MatrixXd &zijOut, Eigen::MatrixXd &lamOut){
    
}

void hardycpp::stresskinetic(const Eigen::MatrixXd &InsidersIn, double avvxIn, double avvyIn, double avvzIn, double volIn, Eigen::MatrixXd &SkOut){
    
}

void hardycpp::stresspotential(const Eigen::MatrixXd &FxIn, const Eigen::MatrixXd &FyIn, const Eigen::MatrixXd &FzIn,
                     const Eigen::MatrixXd &xijIn, const Eigen::MatrixXd &yijIn, const Eigen::MatrixXd &zijIn, const Eigen::MatrixXd &lamIn, double volIn, Eigen::MatrixXd &SpOut){
    
}


void hardycpp::test(){
    
    cout<<"Natoms: "<<sdata.natoms<<"\nxmin xmax: "<<sdata.xmin<<" "<<sdata.xmax<<"\nymin ymax: "<<sdata.ymin<<" "<<sdata.ymax<<"\nzmin zmax: "<<sdata.zmin<<" "<<sdata.zmax<<endl;
    
    vector<int> indxs=findindxs(false, 1000, -1.0, 66.0, -1.0, 66.0, 0.0, 0.0);
    
    double *xData=(double*)malloc(sizeof(double)*indxs.size());
    double *yData=(double*)malloc(sizeof(double)*indxs.size());
    for (int i=0; i<indxs.size(); i++) {
        xData[i]=wdata[indxs[i]][7];
        yData[i]=wdata[indxs[i]][8];
        cout<<"Number i="<<i<<" -> "<<indxs[i]<<endl;
    }
    plot(xData, yData, (int)indxs.size());
    delete []yData;
    delete []xData;
}

void hardycpp::plot(const double *xData,const double *yData,int dataSize){
    FILE *gnuplotPipe,*tempDataFile;
    string tempDataFileName;
    double x,y,r;
    tempDataFileName = "tempData";
    gnuplotPipe = popen("/opt/local/bin/gnuplot","w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title \"Flow Region\" \n");
        fprintf(gnuplotPipe,"plot \"%s\" with points \n",tempDataFileName.c_str());
        fflush(gnuplotPipe);
        tempDataFile = fopen(tempDataFileName.c_str(),"w");
        for (int i=0; i <= dataSize; i++) {
            x = xData[i];
            y = yData[i];
            r = 0.3;
            fprintf(tempDataFile,"%lf %lf %lf\n",x,y,r);
        }
        fclose(tempDataFile);
        printf("press enter to continue...");
        getchar();
        remove(tempDataFileName.c_str());
        fprintf(gnuplotPipe,"exit \n");
    } else {
        printf("gnuplot not found...");
    }
}
