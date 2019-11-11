#include "write.h"

_Bool inThere(int some, int* There, int n){
    int i; _Bool flag = 0;
    for(i = 0; i < n; ++i){
        if(some == There[i]){flag = 1;}
    }

    return flag;
}

_Bool isEqueal(int ad, int* p, int* compare, int n){
    _Bool flag = 1;

    flag = flag && inThere(ad, compare, n);
    for(int i = 0; i < n - 1; ++i){
        flag = flag && inThere(p[i], compare, n);
    }

    return flag;
}


_Bool canPar(int learn_node, int* p, int n, int ad, constran* s){
    if(ad == learn_node){
        return 0;
    }

    _Bool flag = !inThere(ad, p, n);

    if(s == NULL){
        return flag;
    }

    while(s){
        if(isEqueal(ad, p, s->p, s->np)){
            flag = 0;
        }
        s = s->next;
    }

    return flag;
}

int* bePar(int* p, int n, int ad){
    int* pp = (int*)malloc(n * sizeof(int));
    if(p == NULL){
        pp[0] = ad;
    }else{
        for(int i = 0; i < n - 1; ++i){
            pp[i] = p[i];
        }
        pp[n - 1] = ad;
    }

    return pp;
}

constran* constranFree(constran* cur){
    if(cur == NULL) return cur;
    
    constran* H = cur;
    constran* now = cur->next;

    while(now){
        if(!(now->have)){
            cur->next = now->next;
            free(now->p);
            free(now);
        }else{
            cur = cur->next;
        }
        now = cur->next;
    }

    if(!(H->have))
        H = H->next;

    return H;
}


void showConstran(constran* s){
    printf("I am in showConstran\n");
    while(s){
        printf("%d %f %f\n", s->np, s->score, s->pen);
        s = s->next;
    }
    printf("showConstran done\n");
}
