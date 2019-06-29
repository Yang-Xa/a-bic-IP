#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int count = 0;
//�Ӻ���ǰѡȡ��ѡ��λ��i������ǰi-1������ѡȡm-1����
//�� 1 2 3 4 5 ��ѡȡ 3 ��
//1��ѡȡ5������ǰ4������ѡȡ2������ǰ4������ѡȡ2������һ�������⣬�ݹ鼴�ɡ�
//2�����������5��ֱ��ѡ��4����ô����ǰ3������ѡȡ2������ǰ��������ѡȡ2������һ�������⣬�ݹ鼴�ɡ�
//3�����Ҳ������4��ֱ��ѡȡ3����ô����ǰ2������ѡȡ2�����պ�ֻ��������
//���򿴣�1��2��3�պ���һ��forѭ������ֵΪ5����ֵΪm
//���򿴣�������Ϊһ��ǰi-1����ѡm-1�ĵݹ顣
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

	for (int i = nLen; i >= m; --i)	//�Ӻ���ǰ����ѡ��һ��
	{
		out[m-1] = arr[i-1]; //ѡ��һ����
		Combination(arr, i-1, m-1, out, outLen); // ��ǰi-1������ѡȡm-1�����еݹ�
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
