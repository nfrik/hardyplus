//
//  main.c
//  lammps2txt
//
//  Created by Nikolay Frik on 6/14/13.
//  Copyright (c) 2013 Nikolay Frik. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>


int main(int argc, const char * argv[])
{    
    char cwf[1024]; // current working file
    char fileNameRead[1024];
    char fileNameWrite[1024];
    
    
    if(getcwd(cwf, sizeof(cwf)) == NULL){
        perror("getcwd(...) error");
    }
    
    if(argc<2){
        printf("%s","Error: missing path to file \n ");
        return 0;
    }
    
    if(argv[1][0]=='.')
        strcat(cwf, &argv[1][1]);
    else if(argv[1][0]=='/') //we assume user supplies full path
    {
        strcpy(cwf, &argv[1][0]);
    }
    else //we assume user
    {
        strcat(cwf,"/");
        strcat(cwf, &argv[1][0]);
    }
        
    printf("Current working file: %s \n ",cwf);
    
    //    int strl = (int)strlen(argv[1]);
    //    char *fileNameRead = malloc(sizeof(char)*strl);
    //    char *fileNameWrite = malloc(sizeof(char)*strl);
    
    strcpy(fileNameRead, cwf); //copy read filename first
    strcpy(fileNameWrite, fileNameRead); //copy write filename
    
    int strl = (int)strlen(fileNameRead);

    fileNameWrite[strl-4]='.';
    fileNameWrite[strl-3]='x';
    fileNameWrite[strl-2]='t';
    fileNameWrite[strl-1]='x';
    
    //main function changes extension letters only
    
    printf("Input file: %s \n ",fileNameRead);
    printf("Output file %s \n ",fileNameWrite);
    
    FILE *fileRead = fopen(fileNameRead, "r");
    FILE *fileWrite = fopen(fileNameWrite, "w");
    
    if (fileRead) {
        char line[BUFSIZ];
        char timestepline[16];
        char linebuf[256];
        int z=0;
        /*        for (int i=0; i<15; i++) {
         fgets(line, sizeof(line), file);
         fputs(line, stdout);
         }
         */
        while (fgets(line, sizeof(line), fileRead)) { //read through all lines
            if (strcmp("ITEM: TIMESTEP",line)<=0) { //find in lines pattern ITEM: TIMESTEP
                fgets(timestepline, sizeof(line), fileRead); //read next line to get timestep
                for (int k=1; k<8+z; k++) {//skip next 8 lines
                    fgets(line, sizeof(line), fileRead);
                    if ((z<1)&&(k<7)) {
                        fputs(line, fileWrite);
                    }
                }
            }//if
            
            
            //now need to remove \n character from the first string
            char colon = '\n';
            char *found=strchr(line, colon);
            if (found) {
                size_t len = found - line;
                line[len]='\0';
            }
            
            if (z<1){
                snprintf(linebuf, sizeof linebuf, "%s%s",line,"timestep\n");
                z=1;
            }
            else
                snprintf(linebuf, sizeof linebuf, "%s%s",line,timestepline);
            
            fputs(linebuf, fileWrite);
            // fgets(line, sizeof(line), fileRead);//read next line
            
        }
        fclose(fileRead);
        fclose(fileWrite);
    }
    else{
        perror(fileNameRead);
        perror(fileNameWrite);
    }
    
    //free(fileNameRead);
    //free(fileNameWrite);
    
    // insert code here...
    //    printf("Hello, World!\n");
    return 0;
}

