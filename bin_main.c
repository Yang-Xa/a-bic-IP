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
    Parent p; int i,j,k,nm; double base, temp;
    int Pcount;

    time_t start, end;
    time(&start);

    SaveVariables(data.nVariables, argv[3]);
    Save* H;
    Save* cur;
    for(learn_node = 0; learn_node < data.nVariables; ++learn_node){

        int* arr = (int*)malloc((data.nVariables - 1) * sizeof(int));
        for(i = 0, j = 0; i < data.nVariables; ++i){
            if(i == learn_node){continue;}
            arr[j++] = i;
        }

        H = (Save*)malloc(sizeof(Save));
        cur = H;
        //¸¸¼¯ÎªNULL
        p.nParent = 0;
        p.parents = NULL;
        base = score(learn_node, &p, alpha, &data);
        //SaveScore(learn_node, p, base, argv[3]);
        (cur->next) = addSave(base, p);
        cur = cur->next;

        Pcount = 1;

        for(nm = 1; nm <= m; ++nm){
            int* pTable = (int*)malloc((data.nVariables - 1) * sizeof(int));
            memset(pTable, 0, (data.nVariables - 1) * sizeof(int));

            for(i = 0; i < nm; ++i){
                pTable[i] = 1;
            }

            int bFind = 0;
            do{
                p.nParent = nm;
                p.parents = (int*)malloc(p.nParent * sizeof(int));
                for(i = 0, j = 0; i < data.nVariables - 1; ++i){
                    if(pTable[i]){(p.parents)[j++] = arr[i];}
                }
                temp = score(learn_node, &p, alpha, &data);
                if(temp > base){
                    cur->next = addSave(temp, p);
                    cur = cur -> next;

                    ++Pcount;
                }
                free(p.parents);

                bFind = 0;
                for(i = 0; i < data.nVariables - 2; ++i){
                    if(pTable[i] == 1 && pTable[i + 1] == 0){
                        swap(pTable, i, i + 1);
                        bFind = 1;

                        if(!pTable[0]){
                            for(k = 0, j = 0; k < i; ++k){
                                if(pTable[k]){swap(pTable, k, j); ++j;}
                            }
                        }
                        break;
                    }
                }
            }while(bFind);
            free(pTable);
        }


        SaveNode(H, learn_node, Pcount, argv[3]);
        free(H);free(arr);
        printf("have learned node %d\n", learn_node);
    }

    time(&end);
    printf("spend %d s\n", end-start);

    return 0;
}
