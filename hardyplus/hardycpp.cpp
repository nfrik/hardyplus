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
#include <Eigen/Core>

hardycpp::hardycpp(const char *path){
    
    
    this->readRawDataFromFile(path);//read file to wdata
    
    this->run(10000,10,10,1);//run 10000'th frame with 10x10x1 mesh grid
    
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
    
    vector<int> headIndxs=findindxs(true, time, 1.0-rcx, 1.0, -INFINITY, INFINITY, -INFINITY, INFINITY);
    vector<int> tailIndxs=findindxs(true, time, 0.0, rcx, -INFINITY, INFINITY, -INFINITY, INFINITY);
    
    
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

void hardycpp::plot(double *xData,double *yData,int dataSize){
    FILE *gnuplotPipe,*tempDataFile;
    string tempDataFileName;
    double x,y,r;
    int i;
    tempDataFileName = "tempData";
    gnuplotPipe = popen("/opt/local/bin/gnuplot","w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title \"Flow Region\" \n");
        fprintf(gnuplotPipe,"plot \"%s\" with circles lw 0.3\n",tempDataFileName.c_str());
        fflush(gnuplotPipe);
        tempDataFile = fopen(tempDataFileName.c_str(),"w");
        for (i=0; i <= dataSize; i++) {
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
