#include "score.h"
#include <string.h>

typedef struct save{
    int n;
    int* p;
    double score;
    struct save* next;
}Save;

typedef struct constran{
    int np;
    int* p;
    double score;
    double pen;
    _Bool ingore;
    _Bool have;
    struct constran* next;
}constran;


void swap(int* arr, int a, int b){
    int c = *(arr + a);
    *(arr + a) = *(arr + b);
    *(arr + b) = c;
}


Save* addSave(double score, Parent p){
    Save* s = (Save*)malloc(sizeof(Save));
    (s->n) = p.nParent;
    //(s->p) = p.parents;
    (s->p) = (int*)malloc((s->n)*sizeof(int));
    memcpy(s->p, p.parents, (s->n)*sizeof(int));
    (s->score) = score;
    (s->next) = NULL;

    return s;
}

constran* addConstran(Parent p, double score, double pen){
    constran* s = (constran*)malloc(sizeof(constran));
    s->np = p.nParent;
    s->p = (int*)malloc(p.nParent * sizeof(int));
    memcpy(s->p, p.parents, p.nParent * sizeof(int));
    s->score = score;
    s->pen = pen;
    s->have = 1;
    s->ingore = 0;
    s->next = NULL;

    return s;
}

int writeConstran(int learn_node, constran* Group[], int n, char* where){
    int i, j; constran* cur, *tmp; int com = 0;
    for(i = 0; i <= n; ++i){
        cur = Group[i];
        while(cur){
            ++com;
            cur = cur->next;
        }
    }
    // printf("%d have %d parent set\n", learn_node, com);
    FILE *wf = fopen(where, "a+");
    if(wf == NULL){printf("can not open the file\n"); return 0;}
    fprintf(wf, "%d %d\n", learn_node, com);

    for(i = 0; i <= n; ++i){
        cur = Group[i];
        while(cur){
            fprintf(wf, "%f %d ", cur->score, cur->np);
            for(j = 0; j < cur->np; ++j)
                fprintf(wf, "%d ", (cur->p)[j]);
            fprintf(wf, "\n");

            tmp = cur;
            cur = cur->next;
            if(!tmp){free(cur->p);}
            free(tmp);
        }
    }

}

int SaveNode(Save* s, int learn_node, int Pcount, char* where){
    FILE *wf = fopen(where, "a+");
    if(wf == NULL){printf("can no open write file\n"); return 0;}
    int i;
    fprintf(wf, "%d %d\n", learn_node, Pcount);

    Save* cur;
    s = s->next;
    while(s != NULL){
        cur = s;
        s = s->next;

        fprintf(wf, "%f %d ", cur->score, cur->n);
        for(i = 0; i < cur->n; ++i)
            fprintf(wf, "%d ", (cur->p)[i]);
        fprintf(wf, "\n");

        if((cur->p) != NULL){free(cur->p);}
        free(cur);
    }

    fclose(wf);
    return 1;
}


int SaveVariables(int n, char* where){
    FILE *wf = fopen(where, "a+");
    if(wf == NULL){printf("can no open write file\n"); return 0;}

    fprintf(wf, "%d\n", n);

    fclose(wf);
    return 1;
}


/*
int SaveScore(int node, Parent p, double score, char* where){
    FILE *wf = fopen(where, "a+");
    if(wf == NULL){printf("can no open write file\n"); return 0;}

    //if(!(p.nParent)){
    //    fprintf(wf, "%d\n", node);
    //}

    //fprintf(wf, "%d\t%f\t%d\t", node, score, p.nParent);
    fprintf(wf, "%f %d ", score, p.nParent);
    int i;
    for(i = 0; i < p.nParent; ++i)
        fprintf(wf, "%d ", p.parents[i]);
    fprintf(wf, "\n");

    fclose(wf);
    return 1;
}
*/
