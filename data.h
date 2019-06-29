#include <stdlib.h>
#include <time.h>

typedef struct Data{
    int nSamples;
    int nVariables;
    char* data;
    int* arities;
}Data;

void getArity(Data* mdata){
    (mdata -> arities) = (int*)malloc((mdata -> nVariables) * sizeof(int));
    int i, j;
    for(i = 0; i < (mdata -> nVariables); ++i){
        int val[3] = {0};
        for(j = 0; j < (mdata -> nSamples); ++j){
            if(!(val[*((mdata -> data) + j * (mdata -> nVariables) + i) - '0'] ))
                val[*((mdata -> data) + j * (mdata -> nVariables) + i) - '0'] = 1;
        }
        *((mdata -> arities) + i) = val[0] + val[1] + val[2];
    }
}

int getData(char* fname, Data* mdata){
    FILE* fp = fopen(fname, "r");
    if(fp == NULL){
        printf("can't open the file\n");
        return 0;
    }
    int c;
    (mdata -> nVariables) = 1;
    while((c = getc(fp)) != '\n'){
        if(c == '\t')
            ++(mdata -> nVariables);
    }
    (mdata -> nSamples) = 0;
    while((c = getc(fp)) != EOF){
        if(c == '\n'){++(mdata -> nSamples);}
    }

    fseek(fp, 0,SEEK_SET);
    while((c = getc(fp)) != '\n'){}
    (mdata -> data) = (char*)malloc((mdata -> nSamples) * (mdata -> nVariables) * sizeof(char));
    int n = 0;
    while((c = getc(fp)) != EOF){
        if(c != '\t' && c != '\n' && c != '\r'){
            *((mdata -> data) + n) = c;
            ++n;
        }
    }
    fclose(fp);
    getArity(mdata);
    return 1;
}


/*Ð§ÂÊÂý
    (data ->nSamples) = 1;
    (data -> data) = (char*)malloc((data -> nVariables) * sizeof(char));
    while((c = getc(fp)) != EOF){
        int n = 0;
        while(n < (data -> nVariables)){
            if(c != '\t' && c != '\n'){
                *((data -> data) + n * ((data->nSamples)-1)) = c;
                ++n;
            }
        }
        ++(data -> nSamples);
        (data -> data) = (char*)realloc((data -> data), (data -> nSamples)*(data -> nVariables)*sizeof(char));
    }
*/
