#include <stdio.h>
#include "write.h"


int* Combination(int* arr, int nLen, int m, int* out, int outLen)
{
	if(m == 0)
	{
		return out;
	}

	for (int i = nLen; i >= m; --i)	//从后往前依次选定一个
	{
		out[m-1] = arr[i-1]; //选定一个后
		Combination(arr, i-1, m-1, out, outLen); // 从前i-1个里面选取m-1个进行递归
	}
}

int main(int argc, char* argv[]){
    if(!argv[1]){
        printf("usage: ./main inputfile\n");
        return -1;
    }

    Data data;
    if(!getData(argv[1], &data)){return 0;}

    //int learn_node = atoi(argv[2]);
    int m = atoi(argv[4]);
    double alpha = atof(argv[2]);
    int learn_node;
    Parent p; int i,j,k; double base, temp;

    time_t start, end;
    time(&start);

    for(learn_node = 0; learn_node < data.nVariables; ++learn_node){

        //父集为NULL
        p.nParent = 0;
        p.parents = NULL;
        base = score(learn_node, &p, alpha, &data);
        SaveScore(learn_node, p, base, argv[3]);

        i = 1;
        while(i <= m){
            p.nParent = i;
            int* arr = (int*)malloc((data.nVariables - 1) * sizeof(int));
            for(j = 0, k = 0; j < data.nVariables; ++j){
                if(j == learn_node){continue;}
                arr[k++] = j;
            }
            int* out = (int*)malloc(i * sizeof(int));
            p.parents = Combination(arr, nLen, i, out, i);

            free(arr);free(out);
            ++i;
        }



    }

    time(&end);
    printf("spend %d s\n", end-start);

    return 0;
}

/*

//父集个数为1
        if(m >= 1){
            p.nParent = 1;

            p.parents = (int*)malloc(p.nParent * sizeof(int));
            for(i = 0; i < data.nVariables; ++i){
                if(i == learn_node){continue;}
                p.parents[0] = i;
                temp = score(learn_node, &p, alpha, &data);
                if(temp > base){SaveScore(learn_node, p, temp, argv[3]);}
            }
            free(p.parents);
        }
*/
