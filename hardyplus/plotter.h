//
//  plotter.h
//  hardyplus
//
//  Created by Nikolay Frik on 9/11/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#ifndef __hardyplus__plotter__
#define __hardyplus__plotter__

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#endif /* defined(__hardyplus__plotter__) */


using namespace std;

class plotter{
    public:
    plotter();
    ~plotter();
    void open();
    void plot(double *xData,double *yData,int dataSize);
    void clear();
    void close();
    void gen_random(char *s, const int len);
    
    private:
    FILE *gnuplotPipe,*tempDataFile;
    string tempDataFileName;
};
