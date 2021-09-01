#include <stdio.h>
#include <stdlib.h>
int **ano_approach(int **h, int n, int m);

int **ano_approach(int **h, int n, int m){
    int **vec1=(int **)malloc(sizeof(int*)*(n+1));
    int **vec2=(int **)malloc(sizeof(int*)*(m+1));
    int **vec3=(int **)malloc(sizeof(int*)*(m+1));
    for(int i=0;i<n+1;i++) vec1[i]=(int *)malloc(sizeof(int)*(m+1));
    for(int i=0;i<m+1;i++){
        vec2[i]=(int *)malloc(sizeof(int)*(m+1));
        vec3[i]=(int *)malloc(sizeof(int)*(m+1));
    }

    for(int i=1;i<n+1;i++){
        memset(vec1[i],0,sizeof(int)*(m+1));
        int count1=0;
        for(int j=1;j<m+1;j++){
            if(h[i][j]==1)  vec1[i][++count1]=j;
        }
    }

    for(int i=1;i<m+1;i++) memset(vec2[i],0,sizeof(int)*(m+1));
    for(int i=1;i<m+1;i++) vec2[i][i]=1;
    for(int i=1;i<n+1;i++){
        int temp2=vec1[i][1];
        vec2[temp2][temp2]=0;
        int count2=2;
        while(vec1[i][count2]>0&&count2<m+1){
            int a=vec1[i][count2];
            vec2[temp2][a]=1;
            count2++;
        }
    }

    int count3=1;
    for(int j=1;j<m+1;j++){
        int temp3=0;
        for(int i=1;i<m+1;i++)
            temp3+=vec2[i][j];
        if(temp3==0) continue;
        else{
            for(int k=1;k<m+1;k++)
                vec3[count3][k]=vec2[k][j];
            count3++;
        }
    }

    for(int i=0;i<n+1;i++)  free(vec1[i]);
    for(int i=0;i<m+1;i++)  free(vec2[i]);
    free(vec1);free(vec2);
    return vec3;
}
