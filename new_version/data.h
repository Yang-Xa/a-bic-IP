#include <stdlib.h>
#include <time.h>
#include <string.h>

typedef struct Data{
    int nSamples;
    int nVariables;
    char* data;
    int* arities;
    char** varNames;
}Data;

//candidate parent set
typedef struct parCan{
    int ** parentCan;
}parCan;

void someFree(Data mdata, parCan par){
    if(mdata.data != NULL) free(mdata.data);
    if(mdata.arities != NULL) free(mdata.arities);

    int i;
    for(i = 0; i < mdata.nVariables; ++i){
        if(par.parentCan[i] != NULL){
            free(par.parentCan[i]);
        }
    }
    free(par.parentCan);
}

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
        printf("%s can not open\n", fname);
        return 0;
    }
    int c;
    (mdata -> nVariables) = 1;
    while((c = getc(fp)) != '\n'){
        if(c == '\t') ++(mdata -> nVariables);
    }
    (mdata -> nSamples) = 0;
    while((c = getc(fp)) != EOF){
        if(c == '\n'){++(mdata -> nSamples);}
    }

    fseek(fp, 0, SEEK_SET);
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

int getParCan(char *fname, parCan *par, Data data){
    (par->parentCan) = (int **)malloc(data.nVariables * sizeof(int *));
    int i;
    for(i = 0; i < data.nVariables; ++i){
        (par->parentCan)[i] = NULL;
    }

    FILE *fp = fopen(fname, "r");
    if(fp == NULL){
        printf("%s can not open\n", fname);
        return 0;
    }

    char *tmp = (char *)malloc(100000 * sizeof(char));
    memset(tmp, 0, 100000);

    int ff;
    int node, tmpnode;

    while(fgets(tmp, 100000, fp) != NULL){
        ff = 0;
        for(i = 0; tmp[i] != '-'; ++i){
            if(tmp[i] != ' '){
                tmpnode = 0;
                while(tmp[i] != ' '){
                    tmpnode = tmpnode * 10 + (tmp[i] - '0');
                    ++i;
                }
                if(ff == 0){node = tmpnode;}
                else if(ff == 1){(par->parentCan)[node] = (int *)malloc((tmpnode + 1) * sizeof(int)); (par->parentCan)[node][ff - 1] = tmpnode;}
                else{(par->parentCan)[node][ff - 1] = tmpnode;}
                ++ff;
            }
        }
    }

    fclose(fp);

    return 1;
}

/*Ч����
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
