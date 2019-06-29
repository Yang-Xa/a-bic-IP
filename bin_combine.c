#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void swap(int* arr, int a, int b)
{
    int c = *(arr + a);
    *(arr + a) = *(arr + b);
    *(arr + b) = c;
}

int Combine(int* arr, int n, int m)
{
    if(m > n){return 0;}

    int comm = 0;

    int* pTable = (int*)malloc(n * sizeof(int));
    memset(pTable, 0, n * sizeof(int));

    int i, j, k;

    for(i = 0; i < m; ++i)
        pTable[i] = 1;

    int bFind = 0;
    do
    {
        for(i = 0; i < n; ++i)
        {
            if(pTable[i])
                printf("%d ", arr[i]);
        }
        printf("\n");

        ++comm;

        bFind = 0;
        for(i = 0; i < n - 1; ++i)
        {
            if(pTable[i] == 1 && pTable[i + 1] == 0)
            {
                swap(pTable, i, i + 1);
                bFind = 1;

                if(!pTable[0])
                {
                    for(k = 0, j = 0; k < i; ++k)
                    {
                        if(pTable[k]){swap(pTable, k, j); ++j;}
                    }
                }
                break;
            }
        }
    }while(bFind);

    free(pTable);
    return comm;
}


int main(int argc, char* argv[])
{
    int arr[6] = {1, 2, 3, 4, 5, 6};

    int aaa[1000];
    for(int i = 0; i < 1000; ++i)
        aaa[i] = i + 1;

    time_t start, end;

    time(&start);
    printf("have %d things\n", Combine(aaa, 1000, 2));
    time(&end);

    printf("spend %d s\n", end - start);

    return 0;
}
