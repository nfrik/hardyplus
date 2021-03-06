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
//#include "plotter.h"

using namespace Eigen;

template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

hardycpp::hardycpp(const char *path){
    
    this->sdata.timestep=1000;
    this->sdata.nframes=0;
    this->workFilePath=(char *)path;
    //this->readRawDataFromFile(path);//read file to wdata
    this->openWorkFile();

    //we specify current folder path by truncating file name
    folderpath=path;
    unsigned long last=folderpath.find_last_of("/\\");
    folderpath=folderpath.substr(0,last+1);
}

hardycpp::~hardycpp(){
    this->closeWorkFile();
}

/**
 * Calculate P12 for one frame
 *
 */
void hardycpp::run(int time, int dargx, int dargy, int dargz){
//    cout<<"Natoms: "<<sdata.natoms<<"\nxmin xmax: "<<sdata.xmin<<" "<<sdata.xmax<<"\nymin ymax: "<<sdata.ymin<<" "<<sdata.ymax<<"\nzmin zmax: "<<sdata.zmin<<" "<<sdata.zmax<<endl;
    
    //wee need to compensate atom coordinates because of removal of walls
//    double yclo=0.06675+0.00001;
//    double ychi=0.9166-0.00001;
    double yclo=0.04009+0.00001;
    double ychi=0.9599-0.00001;
    double dely=ychi-yclo;
    
    
    //Rescale rc (potential cutoff) from lj units and true box dimension
    //double rc=1.12246;//potential cutoff distance
    double rc=2.2;//potential cutoff distance
    double rcx=rc/(sdata.xmax-sdata.xmin);
    double rcy=rc/(sdata.ymax-sdata.ymin);
    double rcz=rc/(sdata.zmax-sdata.zmin);
    
    double sxlo,sylo,szlo,sxhi,syhi,szhi,dx,dy,dz,xlo,xhi,ylo,yhi,zlo,zhi,dsx,dsy,dsz;
    double tmass,tvx,tvy,tvz, tpx, tpy, tpz, rho,vol;
    long    tN; //number of particles
    
    //data preprocessing
    //pad data with tail and head in x direction
    //1***************2%     %2***************1%
    //1***************2%     %2***************1%
    //1***************2% --> %2***************1%
    //1***************2%     %2***************1%
    //1***************2%     %2***************1%
    //---------------------------
    
    //plotter plotter;
    
    MatrixXd data, U, Fx, Fy, Fz, xij, yij, zij, lam, inatoms, outatoms, Sk, Sv, output;
    
    
    output.resize(dargy, dargx);
    
    
    if(!this->readNextFrame()){
        cout << "Something went wrong with work file in hardy::run(...) ";
    }
    
    data.setConstant(0);
    //extract from wdata[][] -> data
    getBodyHeadTail2Matrix(data, rcx);
    
    hardycpp::printMat2File(data, folderpath+"sample_m.txt");
    
    for (int i=1; i<=dargy; i++) {
        for (int j=1; j<=dargx; j++) {
            for (int k=1; k<=dargz; k++) {
                
                U.setConstant(0);
                Fx.setConstant(0); Fy.setConstant(0); Fz.setConstant(0);
                xij.setConstant(0); yij.setConstant(0); zij.setConstant(0);
                lam.setConstant(0);
                inatoms.setConstant(0); outatoms.setConstant(0);
                Sk.setConstant(0); Sv.setConstant(0);
                
                sxlo=(j-1.0)/dargx;
                sxhi=j*1.0/dargx;
                sylo=(i-1.0)*dely/dargy+yclo;
                syhi=i*dely*1.0/dargy+yclo;
                szlo=(k-1.0)/dargz;
                szhi=k*1.0/dargz;
                
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
                
                //vector<int> inatoms=findindxs(true, time, sxlo, sxhi, sylo, syhi, szlo, szhi);
                //vector<int> outatoms=findindxs(true, time, sxlo-rcx, sxhi+rcx, sylo-rcy, syhi+rcy, szlo-rcz, szlo+rcz);

                getInsideAtoms(data, inatoms, true, sxlo, sxhi, sylo, syhi, -INFINITY, INFINITY);
                getOutsideAtoms(data, outatoms, true, sxlo, sxhi, sylo, syhi, -INFINITY, INFINITY, rcx, rcy, rcz);
                
                neighborList(inatoms, outatoms, U, Fx, Fy, Fz, xij, yij, zij, lam,rc);
                
                tN=inatoms.col(10).size();   //total number of particles
                tmass=inatoms.col(10).sum(); //total mass
                tvx=inatoms.col(4).sum()/tN; //average velocity
                tvy=inatoms.col(5).sum()/tN; //average velocity
                tvz=inatoms.col(6).sum()/tN; //average velocity
                rho=tmass/vol; //density
                
                tpx=tmass*tvx/vol;
                tpy=tmass*tvy/vol;
                tpz=tmass*tvz/vol;
                
                stresskinetic(inatoms, tvx, tvy, tvz, vol, Sk);
                stresspotential(Fx, Fy, Fz, xij, yij, zij, lam, vol, Sv);
                
                output(i-1,j-1)=-(Sv(0,1)+Sk(0,1));
                
                printMat2File(U, folderpath+string("Udat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(Fx, folderpath+string("Fxdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(Fy, folderpath+string("Fydat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(Fz, folderpath+string("Fzdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(xij, folderpath+string("Xijdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(yij, folderpath+string("Yijdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(zij, folderpath+string("Zijdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(lam, folderpath+string("Lamdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(inatoms, folderpath+string("InAtomsdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(outatoms, folderpath+string("OutAtomsdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(Sk, folderpath+string("Skdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                printMat2File(Sv, folderpath+string("Svdat_i")+NumberToString(i)+string("_j_")+NumberToString(j)+string(".txt"));
                
                //cout<<"Pxy for cell["<<i<<"]["<<j<<"] = "<<(Sk+Sv).row(0).col(0)<<endl;
                
                
                
                //neighborList(inatoms, outatoms, U, Fx, Fy, Fz, xij, yij, zij, lam);
                
            }
        }
    }
    
    printMat2File(output, folderpath + string("P12_t_")+NumberToString(time)+string(".txt"));

}

void hardycpp::openWorkFile(){
    if(!this->workFileStream.is_open()){
     this->workFileStream.open(this->workFilePath);
    }
    
    string line;
    
    //determine number of lines in file
    int flines=0;
    
    if (this->workFileStream.is_open())
        while (getline(this->workFileStream, line))
            flines++;
    
    this->workFileStream.clear();
    this->workFileStream.seekg(0,ios::beg);
    
    if (this->workFileStream.is_open()) {
        //read preamble
        for (int i=0; i<7; i++){
            getline(this->workFileStream, line);
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
    }
    
    sdata.nframes=flines/sdata.natoms;
    
    //resize our matrix for
    this->wdata.resize(sdata.natoms); // set height
    for(int i=0;i<sdata.natoms;i++)
        this->wdata[i].resize(12); // set width
}

void hardycpp::closeWorkFile(){
    if(this->workFileStream.is_open()){
        this->workFileStream.close();
    }
}

bool hardycpp::readNextFrame(){
    
    string line;
    int i=0;
    for (i=0; i<sdata.natoms; i++) {
        if (!getline(this->workFileStream, line)) {
            cout << "EOF or something is wrong";
            return false;
        }
        sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&wdata[i][0],&wdata[i][1],&wdata[i][2],&wdata[i][3],&wdata[i][4],&wdata[i][5],&wdata[i][6],&wdata[i][7],&wdata[i][8],&wdata[i][9],&wdata[i][10],&wdata[i][11]);
    }
    cout << "Finished reading " << wdata[i-1][11]<<"'th frame";
    return true;
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

void hardycpp::findatoms(Eigen::MatrixXd &atoms,bool scaled, int time, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax){
    vector<int> indexes=findindxs(scaled, time, xmin, xmax, ymin, ymax, zmin, zmax);
    long rows=indexes.size();
    atoms.resize(rows, 11);
    
    for (int i=0; i<rows; i++) {
        atoms.row(i)=VectorXd::Map(&wdata[indexes[i]][0], 11);
    }

 }

//Function assembles particles within cutoff range from left and right borders and returns it through matrix argument
//Matrix returned through input argument m,
//rcx - is a cutoff distance of lj potential
void hardycpp::getBodyHeadTail2Matrix(Eigen::MatrixXd &m,double rcx){
    vector<int> indexHead;
    vector<int> indexTail;
    int i; //timestep by default is scale of 1000
    
//    vector<int> headIndxs=findindxs(true, time, 1.0-rcx, 1.0, -INFINITY, INFINITY, -INFINITY, INFINITY);
//    vector<int> tailIndxs=findindxs(true, time, 0.0, rcx, -INFINITY, INFINITY, -INFINITY, INFINITY);
    
    
    for (i=0; i<sdata.natoms; i++) {
        if ((1.0-rcx<=wdata[i][1])&&(wdata[i][1]<=1.0))
            indexHead.push_back(i);
        else if((0.0<=wdata[i][1])&&(wdata[i][1]<=rcx))
            indexTail.push_back(i);
    }

    
    m.resize(sdata.natoms+indexHead.size()+indexTail.size(), 11);
    
    int j=0;
    
//    head(:,2)=head(:,2)-1;
//    head(:,8)=head(:,8)-(sdata(1,2)-sdata(1,1));
//    tail(:,2)=tail(:,2)+1;
//    tail(:,8)=tail(:,8)+(sdata(1,2)-sdata(1,1));
    
    VectorXd headtailFix(11);
    headtailFix<<0,1,0,0,0,0,0,(sdata.xmax-sdata.xmin),0,0,0;
    
    //piecewize initialization
    for (unsigned long i=0; i<sdata.natoms; i++) {
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
//        cout<<atoms.row(i)<<endl;
    }
    
}

void hardycpp::getOutsideAtoms(const Eigen::MatrixXd &data, Eigen::MatrixXd &atoms, bool scaled, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double rcx, double rcy, double rcz){
    
    vector <int> indexes;
    long end=data.rows();
    int i; //timestep by default is scale of 1000    
    
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
//        cout<<atoms.row(i)<<endl;
    }
    
}

void hardycpp::neighborList(const Eigen::MatrixXd &InsidersIn, const Eigen::MatrixXd &OutsidersIn, Eigen::MatrixXd &phiOut, Eigen::MatrixXd &FxOut, Eigen::MatrixXd &FyOut, Eigen::MatrixXd &FzOut, Eigen::MatrixXd &xijOut, Eigen::MatrixXd &yijOut, Eigen::MatrixXd &zijOut, Eigen::MatrixXd &lamOut, double rc){

    
    int NInsiders=InsidersIn.col(1).size();
    int NOutsiders=OutsidersIn.col(1).size();
    int trows=NInsiders+NOutsiders;

    phiOut.resize(trows,trows);
    FxOut.resize(trows, trows);
    FyOut.resize(trows, trows);
    FzOut.resize(trows, trows);
    xijOut.resize(trows, trows);
    yijOut.resize(trows, trows);
    zijOut.resize(trows, trows);
    lamOut.resize(trows, trows);

    phiOut.setConstant(0);
    FxOut.setConstant(0);
    FyOut.setConstant(0);
    FzOut.setConstant(0);
    xijOut.setConstant(0);
    yijOut.setConstant(0);
    zijOut.setConstant(0);
    lamOut.setConstant(0);
    
//    double rc=2.2;//LJ cutoff potential
    
    Vector3d A(0,0,0);
    Vector3d B(0,0,0);
    Vector3d f(0,0,0);
    Vector3d r(0,0,0);
    
    //calculate interactions within the box
    for (int i=0; i<NInsiders; i++) {
        for (int j=i+1; j<NInsiders; j++) {
            A<<InsidersIn(i,7),InsidersIn(i,8),InsidersIn(i,9);
            B<<InsidersIn(j,7),InsidersIn(j,8),InsidersIn(j,9);
            phiOut(i,j)=interact.LJPotential(A, B, rc, 1.0, 1.0);
            f=interact.LJForce(A, B, rc);
            r=A-B;
            FxOut(i,j)=f(0);
            FyOut(i,j)=f(1);
            FzOut(i,j)=f(2);
            xijOut(i,j)=r(0);
            yijOut(i,j)=r(1);
            zijOut(i,j)=r(2);
            lamOut(i,j)=1;
        }
    }
    
    int m;
    for (int i=0; i<NInsiders; i++) {
        for(int j=0; j<NOutsiders; j++){
            
            m=j+NInsiders;
            A<<InsidersIn(i,7),InsidersIn(i,8),InsidersIn(i,9);
            B<<OutsidersIn(j,7),OutsidersIn(j,8),OutsidersIn(j,9);
            phiOut(i,m)=interact.LJPotential(A, B, rc);
            f=interact.LJForce(A, B, rc); //BIG QUESTION? I need to figure out the order A,B or B,A !!!!!!!!!!!!!!!!
            r=A-B;
            FxOut(i,m)=f(0);
            FyOut(i,m)=f(1);
            FzOut(i,m)=f(2);
            xijOut(i,m)=r(0);
            yijOut(i,m)=r(1);
            zijOut(i,m)=r(2);
            
            if(sqrt(r.array().square().sum())<rc)
                m=m;//calculate lambda
            
            lamOut(i,m)=0;
        }
    }
    
//    u=u+u'; %potential is symmetric
//    fx=fx-fx'; %force is antysymmetric
//    fy=fy-fy'; %force is antysymmetric
//    fz=fz-fz'; %force is antysymmetric
//    rijx=rijx-rijx'; %distance is antysymmetric
//    rijy=rijy-rijy'; %distance is antysymmetric
//    rijz=rijz-rijz'; %distance is antysymmetric
//    lam=lam+lam';
    
     phiOut=(phiOut+phiOut.transpose()).eval();
     FxOut=(FxOut-FxOut.transpose()).eval();
     FyOut=(FyOut-FyOut.transpose()).eval();
     FzOut=(FzOut-FzOut.transpose()).eval();
     xijOut=(xijOut-xijOut.transpose()).eval();
     yijOut=(yijOut-yijOut.transpose()).eval();
     zijOut=(zijOut-zijOut.transpose()).eval();
     lamOut=(lamOut+lamOut.transpose()).eval();
}

double hardycpp::lambda(const Eigen::MatrixXd &iAt,const Eigen::MatrixXd &oAt,const Eigen::MatrixXd &boxDim){
    
    double xlo=boxDim(0,0);
    double xhi=boxDim(0,1);
    double ylo=boxDim(1,0);
    double yhi=boxDim(1,1);
    double mx=0,my=0,mz=0,lam=0,k=0,b=0;
    double ym1=0, xm1=0, ym2=0, xm2=0,ym3=0, xm3=0, ym4=0, xm4=0;
    double eps=numeric_limits<float>::epsilon();
    
    int l=0;
    
    //%check for absolute vertical and horizontal
    if (abs(iAt(0)-oAt(0))<eps) {
        if (oAt(1)>=yhi) {
            mx=iAt(0);
            my=yhi;
            lam=abs(yhi-iAt(1))/abs(iAt(1)-oAt(1));
        }else if(oAt(1)<=ylo){
            mx=iAt(0);
            my=ylo;
            lam=abs(ylo-iAt(1))/abs(iAt(1)-oAt(1));
        }
        l=1;
    }


    if (abs(iAt(1)-oAt(1))<eps) {
        if (oAt(0)<=xhi) {
            mx=xhi;
            my=oAt(1);
            lam=abs(xhi-iAt(0))/abs(iAt(0)-oAt(0));
        }else if(oAt(0)<=xlo){
            mx=xlo;
            my=oAt(1);
            lam=abs(xlo-iAt(0))/abs(iAt(0)-oAt(0));
        }
        l=1;
    }

    if(l<1){
        k=(oAt(1)-iAt(1))/(oAt(0)-iAt(0));
        b=oAt(1)-k*oAt(0);
        
    }
    
    //%find points where wall intersects with A-B line
    ym1=k*xlo+b;
    xm1=xlo;
    ym2=ylo;
    xm2=(ylo-b)/k;
    ym3=yhi;
    xm3=(yhi-b)/k;
    ym4=k*xhi+b;
    xm4=xhi;
    
    
    if(interact.isPMemberOfRectAB(xm1, ym1, iAt(0), iAt(1), oAt(0), oAt(1))){
        return interact.lam(iAt(0), iAt(1), xm1, ym1, oAt(0), oAt(1));
    }else if(interact.isPMemberOfRectAB(xm2, ym2, iAt(0), iAt(1), oAt(0), oAt(1))){
        return interact.lam(iAt(0), iAt(1), xm2, ym2, oAt(0), oAt(1));
    }else if(interact.isPMemberOfRectAB(xm3, ym3, iAt(0), iAt(1), oAt(0), oAt(1))){
        return interact.lam(iAt(0), iAt(1), xm3, ym3, oAt(0), oAt(1));
    }else if(interact.isPMemberOfRectAB(xm4, ym4, iAt(0), iAt(1), oAt(0), oAt(1))){
        return interact.lam(iAt(0), iAt(1), xm4, ym4, oAt(0), oAt(1));
    }else{
        return 1.0;
    }
}

void hardycpp::stresskinetic(const Eigen::MatrixXd &InsidersIn, double avvxIn, double avvyIn, double avvzIn, double volIn, Eigen::MatrixXd &SkOut){
//    function Sk=stresskinetic(psiatoms,vel,vol)
//    Sk=zeros(3,3);
//    v=[psiatoms(:,5) psiatoms(:,6) psiatoms(:,7)];
//    m=psiatoms(:,11);
//    vminusu=[v(:,1)-vel(1) v(:,2)-vel(2) v(:,3)-vel(3)];
//    for i=1:3
//        for j=1:3
//            Sk(i,j)=-sum(m.*vminusu(:,i).*vminusu(:,j));
//    end
//    end
//    Sk=Sk/vol;
//    
//    end
    SkOut.resize(3, 3);
    Vector3d uvel(avvxIn,avvyIn,avvzIn);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            //col 11 - mass
            //col 5 - vx nonscaled
            SkOut(i,j)=-InsidersIn.col(10).cwiseProduct((InsidersIn.col(i+4).array()-uvel(i)).matrix()).cwiseProduct((InsidersIn.col(j+4).array()-uvel(j)).matrix()).sum()/volIn;
        }
    }
}

void hardycpp::stresspotential(const Eigen::MatrixXd &FxIn, const Eigen::MatrixXd &FyIn, const Eigen::MatrixXd &FzIn,
                     const Eigen::MatrixXd &xijIn, const Eigen::MatrixXd &yijIn, const Eigen::MatrixXd &zijIn, const Eigen::MatrixXd &lamIn, double volIn, Eigen::MatrixXd &SpOut){
    
//    Sv(1,1)=-0.5*sum(sum(sum(Fx.*xij.*lam))); %xx
//    Sv(1,2)=-0.5*sum(sum(sum(Fx.*yij.*lam))); %xy
//    Sv(1,3)=-0.5*sum(sum(sum(Fx.*zij.*lam))); %xz
//    
//    Sv(2,1)=-0.5*sum(sum(sum(Fy.*xij.*lam))); %yx
//    Sv(2,2)=-0.5*sum(sum(sum(Fy.*yij.*lam))); %yy
//    Sv(2,3)=-0.5*sum(sum(sum(Fy.*zij.*lam))); %yz
//    
//    Sv(3,1)=-0.5*sum(sum(sum(Fz.*xij.*lam))); %zx
//    Sv(3,2)=-0.5*sum(sum(sum(Fz.*yij.*lam))); %zy
//    Sv(3,3)=-0.5*sum(sum(sum(Fz.*zij.*lam))); %zz
//    
//    Sv=Sv/vol;
    
    int imax=(int)FxIn.rows();
    int jmax=(int)FxIn.cols();
    SpOut.resize(3, 3);
    SpOut.setConstant(0);
    
    for (int i=0; i<imax; i++) {
        for (int j=0; j<jmax; j++) {
            SpOut(0,0)+=FxIn(i,j)*xijIn(i,j)*lamIn(i,j);
            SpOut(0,1)+=FxIn(i,j)*yijIn(i,j)*lamIn(i,j);
            SpOut(0,2)+=FxIn(i,j)*zijIn(i,j)*lamIn(i,j);
            
            SpOut(1,0)+=FyIn(i,j)*xijIn(i,j)*lamIn(i,j);
            SpOut(1,1)+=FyIn(i,j)*yijIn(i,j)*lamIn(i,j);
            SpOut(1,2)+=FyIn(i,j)*zijIn(i,j)*lamIn(i,j);
            
            SpOut(2,0)+=FzIn(i,j)*xijIn(i,j)*lamIn(i,j);
            SpOut(2,1)+=FzIn(i,j)*yijIn(i,j)*lamIn(i,j);
            SpOut(2,2)+=FzIn(i,j)*zijIn(i,j)*lamIn(i,j);
        }
    }
            SpOut=-0.5/volIn*SpOut;
//     SpOut.resize(3, 3);
//     SpOut.setConstant(0);
//     SpOut(0,0)=-0.5*FxIn.cwiseProduct(xijIn).cwiseProduct(lamIn).sum()/volIn;
//     SpOut(0,1)=-0.5*FxIn.cwiseProduct(yijIn).cwiseProduct(lamIn).sum()/volIn;
//     SpOut(0,2)=-0.5*FxIn.cwiseProduct(zijIn).cwiseProduct(lamIn).sum()/volIn;
//
//     SpOut(1,0)=-0.5*FyIn.cwiseProduct(xijIn).cwiseProduct(lamIn).sum()/volIn;
//     SpOut(1,1)=-0.5*FyIn.cwiseProduct(yijIn).cwiseProduct(lamIn).sum()/volIn;
//     SpOut(1,2)=-0.5*FyIn.cwiseProduct(zijIn).cwiseProduct(lamIn).sum()/volIn;
//    
//     SpOut(2,0)=-0.5*FzIn.cwiseProduct(xijIn).cwiseProduct(lamIn).sum()/volIn;
//     SpOut(2,1)=-0.5*FzIn.cwiseProduct(yijIn).cwiseProduct(lamIn).sum()/volIn;
//     SpOut(2,2)=-0.5*FzIn.cwiseProduct(zijIn).cwiseProduct(lamIn).sum()/volIn;
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

void hardycpp::printMat2File(const Eigen::MatrixXd &m, string filename){
    ofstream matoutstr(filename.c_str());
    matoutstr<<m<<endl;
    matoutstr.close();
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
