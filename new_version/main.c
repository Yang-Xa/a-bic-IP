#include <stdio.h>
#include "main.h"

int main(int argc, char* argv[]){
    if(!argv[1] || !argv[2] || !argv[3] || !argv[4] || !argv[5]){
        printf("usage: ./main inputfile_data inputfile_parSet alpha saveWhere maxParent\n");
        return -1;
    }

    // read data
    Data data;
    if(!getData(argv[1], &data)){return 0;}

    // read candidate parent set
    parCan par;
    if(!getParCan(argv[2], &par, data)){return 0;}

    // int learn_node = atoi(argv[2]);
    int m = atoi(argv[5]); int tm;
    //parent set size max is log2N
    int maxm = ceil(log((double) data.nSamples) / log(2.0));
    m = m > maxm ? maxm : m;  

    double alpha = atof(argv[3]);
    int learn_node;
    Parent p; int i, j, k, nm; double base, temp, pen;
    int Pcount;

    time_t start, end;
    time(&start);

    SaveVariables(data.nVariables, argv[4]);


    constran* Group[m + 1];
    constran *preG, *cur;
    int ad;

    // learn all the node's parent
    for(learn_node = 0; learn_node < data.nVariables; ++learn_node){

        // init the Group
        for(i = 0; i < m + 1; ++i)
            Group[i] = NULL;

        // parent set number = 0
        p.nParent = 0;
        p.parents = NULL;
        base = score(learn_node, &p, alpha, &data, &pen);
        Group[0] = addConstran(p, base, pen);

        // 0 < parent set number <= m
        // 候选父集为0
        if(par.parentCan[learn_node] == NULL){
            // 写入父集为0
            writeConstran(learn_node, Group, 0, argv[4]);
            // free(Group[0]);
            printf("have learned node %d\n", learn_node);
            continue;
        }else{
            if(m > par.parentCan[learn_node][0]) tm = par.parentCan[learn_node][0];
            else tm = m;
        }
        // 候选父集不为0
        for(nm = 1; nm <= tm; ++nm){
            preG = Group[nm - 1];
            while(preG){
                if(!(preG->ingore)){
                    for(ad = 1; ad < par.parentCan[learn_node][0] + 1; ++ad){
                        if(canPar(learn_node, preG->p, preG->np, par.parentCan[learn_node][ad], Group[nm])){
                            p.nParent = nm;
                            p.parents = bePar(preG->p, nm, par.parentCan[learn_node][ad]);

                            temp = score(learn_node, &p, alpha, &data, &pen);

                            if(Group[nm] == NULL){
                                Group[nm] = addConstran(p, temp, pen);
                                cur = Group[nm];
                            }else{
                                cur->next = addConstran(p, temp, pen);
                                cur = cur->next;
                            }

                            if(preG->score > temp || temp < base)
                                cur->have = 0;
                            if(data.arities[learn_node] * (cur->pen) + (cur->score) > 0)
                                cur->ingore = 1;
                        }
                    }
                }
                preG = preG->next;
            }


            // showConstran(Group[nm]);
            Group[nm - 1] = constranFree(Group[nm - 1]);
            // printf("%d had learn %d parent set\n", learn_node, nm);
        }

        Group[nm - 1] = constranFree(Group[nm - 1]);

        writeConstran(learn_node, Group, m, argv[4]);

        printf("have learned node %d\n", learn_node);
    }

    time(&end);
    printf("spend %d s\n", end-start);

    someFree(data, par);
    return 0;
}
