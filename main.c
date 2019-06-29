#include <stdio.h>
#include "write.h"

int main(int argc, char* argv[]){
    if(!argv[1] || !argv[2] || !argv[3] || !argv[4]){
        printf("usage: ./main inputfile alpha saveWhere maxParent\n");
        return -1;
    }

    Data data;
    if(!getData(argv[1], &data)){return 0;}

    //int learn_node = atoi(argv[2]);
    int m = atoi(argv[4]);
    double alpha = atof(argv[2]);
    int learn_node;
    Parent p; int i,j,k; double base, temp;
    int Pcount;

    time_t start, end;
    time(&start);

    SaveVariables(data.nVariables, argv[3]);
    Save* H;
    Save* cur;
    for(learn_node = 0; learn_node < data.nVariables; ++learn_node){

        //SaveVariables(0, argv[3]);
        H = (Save*)malloc(sizeof(Save));
        cur = H;
        //父集为NULL
        p.nParent = 0;
        p.parents = NULL;
        base = score(learn_node, &p, alpha, &data);
        //SaveScore(learn_node, p, base, argv[3]);
        (cur->next) = addSave(base, p);
        cur = cur->next;

        Pcount = 1;

        //父集个数为1
        if(m >= 1){
            p.nParent = 1;

            p.parents = (int*)malloc(p.nParent * sizeof(int));
            for(i = 0; i < data.nVariables; ++i){
                if(i == learn_node){continue;}
                p.parents[0] = i;
                temp = score(learn_node, &p, alpha, &data);
                if(temp > base){
                    //SaveScore(learn_node, p, temp, argv[3]);
                    cur->next = addSave(temp, p);
                    cur = cur -> next;

                    ++Pcount;
                }
            }
            free(p.parents);
        }


        //父集个数为2
        if(m >= 2){
            p.nParent = 2;

            p.parents = (int*)malloc(p.nParent * sizeof(int));
            for(i = 0; i < data.nVariables - 1; ++i){
                if(i == learn_node){continue;}
                (p.parents)[0] = i;
                for(j = i + 1; j < data.nVariables; ++j){
                    if(j == learn_node){continue;}
                    (p.parents)[1] = j;
                    temp = score(learn_node, &p, alpha, &data);
                    if(temp > base){
                        //SaveScore(learn_node, p, temp, argv[3]);
                        cur->next = addSave(temp, p);
                        cur = cur->next;

                        ++Pcount;
                    }
                }
            }
            free(p.parents);
        }

        //父集个数为3
        if(m >= 3){
            p.nParent = 3;

            p.parents = (int*)malloc(p.nParent * sizeof(int));
            for(i = 0; i < data.nVariables -2; ++i){
                if(i == learn_node){continue;}
                p.parents[0] = i;
                for(j = i + 1; j < data.nVariables - 1; ++j){
                    if(j == learn_node){continue;}
                    p.parents[1] = j;
                    for(k = j + 1; k < data.nVariables; ++k){
                        if(k == learn_node){continue;}
                        p.parents[2] = k;
                        temp = score(learn_node, &p, alpha, &data);
                        //if(temp > base){SaveScore(learn_node, p, temp, argv[3]); ++Pcount;}
                        if(temp > base){
                            //SaveScore(learn_node, p, temp, argv[3]);
                            cur->next = addSave(temp, p);
                            cur= cur->next;

                            ++Pcount;
                        }
                    }
                }
            }
            free(p.parents);
        }

        SaveNode(H, learn_node, Pcount, argv[3]);
        free(H);

    }

    time(&end);
    printf("spend %d s\n", end-start);


    /*
    Parent p;
    p.nParent = 2;
    int a[2] = {100, 98};
    p.parents = a;

    SaveScore(99, p, score(99, &p, 1.0, &data));
    */

/*
    int ssd[2] = {2, 2};

    printf("%d\n", get3Num(2, ssd));
*/
    /*
    int i;
    for(i = 0; i < data.nVariables; ++i)
        printf("%d  ", (data.arities)[i]);
    */
    //test
    //printf("nVa is %d\n", data.nVariables);
    //printf("nSa is %d\n", data.nSamples);

    return 0;
}
