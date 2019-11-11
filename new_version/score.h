#include "data.h"
#include <math.h>

typedef struct parent{
    int* parents;
    int nParent;
}Parent;

typedef struct combine{
    int combine_jk[3];
}Combine;

//�����ڵ������ת����3����
int get3Num(int n, int* bin){
    int sum = 0;int i;
    for(i = 0; i < n; ++i){
        sum *= 3;
        sum += *(bin + i);
    }
    return sum;
}

double score(int node, Parent* parent, double alpha, Data* mdata, double* mPen){
    double alpha_score = 0.0;

    int r = *((mdata->arities) + node);
    int q = 1, i, j;
    for(i = 0; i < (parent->nParent); ++i)
        q *= *((mdata->arities) + *((parent->parents) + i));

    double Pen = (log((double)(mdata->nSamples))/2)*alpha + 1 - alpha;
    Pen *= (q*(r-1));

    *mPen = Pen;

    int need = (int)pow(3, parent->nParent);
    Combine* com = (Combine*)malloc(need * sizeof(Combine));
    memset(com, 0, need*sizeof(Combine));

    for(i = 0; i < (mdata->nSamples); ++i){
        //int bin[parent->nParent] = {0};
        int* bin = (int*)malloc((parent->nParent) * sizeof(int));
        for(j = 0; j < (parent->nParent); ++j)
            bin[j] = *((mdata->data) + i * (mdata->nVariables) + (parent->parents)[j]) - '0';
        int id = get3Num(parent->nParent, bin);
        ++((com+id)->combine_jk[*((mdata->data) + i * (mdata->nVariables) + node) - '0']);
        free(bin);
    }

    for(i = 0; i < need; ++i){
        int m = (com+i)->combine_jk[0] + (com+i)->combine_jk[1] + (com+i)->combine_jk[2];
        if(!m){continue;}
        for(j = 0; j < 3; ++j){
            if(!(com+i)->combine_jk[j]){continue;}
            alpha_score  += log((double)(com+i)->combine_jk[j]/m) * (double)(com+i)->combine_jk[j];
        }
    }

    free(com);
    return (alpha_score - Pen);
}
