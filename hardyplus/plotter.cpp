//
//  plotter.cpp
//  hardyplus
//
//  Created by Nikolay Frik on 9/11/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#include "plotter.h"
#include <unistd.h>


plotter::plotter(){
    
}

plotter::~plotter(){
    
}

void plotter::open(){
    char rstr[3];
    gen_random(rstr, 3);

    tempDataFileName = string(getcwd(NULL, 0)) + "/tempData" + string(rstr);
    gnuplotPipe = popen("/opt/local/bin/gnuplot","w");
}

void plotter::plot(double *xData,double *yData,int dataSize){
    
//        FILE *gnuplotPipe,*tempDataFile;
//        string tempDataFileName;
        double x,y;
//        tempDataFileName = "tempData";
//        gnuplotPipe = popen("/opt/local/bin/gnuplot","w");
        if (gnuplotPipe) {
            fprintf(gnuplotPipe, "set title \"Flow Region\" \n");
            fprintf(gnuplotPipe,"plot \"%s\" with points \n",tempDataFileName.c_str());
            fflush(gnuplotPipe);
            tempDataFile = fopen(tempDataFileName.c_str(),"w");
            for (int i=0; i <= dataSize; i++) {
                x = xData[i];
                y = yData[i];
                fprintf(tempDataFile,"%lf %lf\n",x,y);
            }
            fclose(tempDataFile);
            //printf("press enter to continue...");
            //getchar();
            usleep(500000);
            remove(tempDataFileName.c_str());

        } else {
            printf("gnuplot not found...");
        }
}

void plotter::clear(){
    if (gnuplotPipe) {        
        fprintf(gnuplotPipe,"clear \n");
    }
    else
        printf("gnuplot not found...");
}

void plotter::close(){
    
    if (gnuplotPipe) {
        fprintf(gnuplotPipe,"exit \n");
    }
    else
        printf("gnuplot not found...");
    
}

void plotter::gen_random(char *s, const int len) {
    srand((unsigned int)time(NULL));    
    static const char alphanum[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
    
    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }
    
    s[len] = 0;
}