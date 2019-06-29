#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int count = 0;
//从后往前选取，选定位置i后，再在前i-1个里面选取m-1个。
//如 1 2 3 4 5 中选取 3 个
//1、选取5后，再在前4个里面选取2个，而前4个里面选取2个又是一个子问题，递归即可。
//2、如果不包含5，直接选定4，那么再在前3个里面选取2个，而前三个里面选取2个又是一个子问题，递归即可。
//3、如果也不包含4，直接选取3，那么再在前2个里面选取2个，刚好只有两个。
//纵向看，1、2、3刚好是一个for循环，初值为5，终值为m
//横向看，该问题为一个前i-1个中选m-1的递归。
void Combination(int arr[], int nLen, int m, int out[], int outLen)
{
	if(m == 0)
	{
		/*for (int j = 0; j < outLen; j++)
			printf("%d\t", out[j]);
		printf("\n");
        */
        ++count;
		return;
	}

	for (int i = nLen; i >= m; --i)	//从后往前依次选定一个
	{
		out[m-1] = arr[i-1]; //选定一个后
		Combination(arr, i-1, m-1, out, outLen); // 从前i-1个里面选取m-1个进行递归
	}
}

void PrintCombination(int arr[], int nLen, int m)
{
	if(m > nLen)
		return;

	int* out = (int*)malloc(m * sizeof(int));
	Combination(arr, nLen, m, out, m);
	free(out);
}


int main(){
    time_t start, end;
    int a[100]; int i;
    for(i = 0; i < 100; ++i)
        a[i] = i;
    time(&start);
    PrintCombination(a, 100, 8);
    time(&end);
    printf("spend %d s\n", end-start);
    printf("count is %d\n", count);

    return 0;
}
