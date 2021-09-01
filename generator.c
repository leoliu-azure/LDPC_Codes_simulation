#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void generate_information_bits(int wanted_length,int *bit_seed, int *bits_process);
void reduced_row_echelon_form(int **row_echelon,int row_count,int column_count);
void add_or_sub(int **row_echelon,int dest, int source, int scalar, int m);
void row_swap(int **row_echelon,int i,int r,int m);
int bin_add(int a, int b);
int bin_multiply(int a, int b);
int **ano_approach(int **h, int n, int m);

int bin_add(int a, int b){
    if((a==1&&b==1)||(a==0&&b==0)) return 0;
    if((a==0&&b==1)||(a==1&&b==0)) return 1;
    return -1;
}

int bin_multiply(int a, int b){
    if(a==1&&b==1) return 1;
    return 0;
}

void generate_information_bits(int wanted_length,int *bit_seed, int *bits_process){
    for(int i=1;i<=6;i++){
        bits_process[i]=bit_seed[i-1];
    }
    for(int i=1;i<=wanted_length-6;i++){
        bits_process[i+6]=bin_add(bits_process[i],bits_process[i+1]);
    }
    int remainder=wanted_length%6;
    for(int i=1;i<=6-remainder;i++){
        bit_seed[i-1]=bits_process[wanted_length+remainder-6+i];
    }
    for(int i=6-remainder+1;i<=6;i++){
        bit_seed[i-1]=bin_add(bits_process[wanted_length+i-6],bits_process[wanted_length+i-5]);
    }
    return;
}

void reduced_row_echelon_form(int **row_echelon,int row_count,int column_count){
    int lead=1;
    column_count+=1;row_count+=1;

   for(int r=1;r<row_count;r++){
        if(column_count<=lead)  return ;
        int i = r;
        while(row_echelon[i][lead]==0){
            i=i+1;
            if(row_count==i){
                i=r;
                lead=lead+1;
                if(column_count==lead)  return ;
            }
        }
        row_swap(row_echelon,i,r,column_count);
        for(i=1;i<row_count;i++){
            if(i!=r){
                int lv=row_echelon[i][lead];
                add_or_sub(row_echelon,i,r,lv,column_count);
            }
        }
        lead=lead+1;
        //printf("r %d lead %d\n",r,lead);
    }
    return ;
}

void add_or_sub(int** row_echelon,int dest, int source, int scalar, int m){
    for(int i=1;i<m;i++){
        row_echelon[dest][i]^=row_echelon[source][i]*scalar;
    }
    return;
}

void row_swap(int** row_echelon,int i,int r,int m){
    if(i==r) return;
    for(int j=1;j<m;j++){
        int temp=row_echelon[i][j];
        row_echelon[i][j]=row_echelon[r][j];
        row_echelon[r][j]=temp;
    }
    return;
}

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

