#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int **identity(int circulant_size);
void right_shift(int** matrix, int circulant_size, int degree, int initial_degree);
int main(const int argc, const char* argv[])
{
    FILE *fin=fopen(argv[1],"r");
    FILE *fout=fopen(argv[2],"w+");

    if(fin==NULL){
        perror("error when opening the file.\n");
        return(-1);
    }
    int bn,bm;
    fscanf(fin,"%d %d",&bn,&bm);
    int** b=(int**)malloc(sizeof(int*)*(bn+1));
    int** mask=(int**)malloc(sizeof(int*)*(bn+1));
    for(int i=0;i<bn+1;i++){
        b[i]=(int*)malloc(sizeof(int)*(bm+1));
        mask[i]=(int*)malloc(sizeof(int)*(bm+1));
    }
    for(int i=1;i<bn+1;i++){
        for(int j=1;j<bm+1;j++){
            fscanf(fin,"%d",&b[i][j]);
        }
    }
    int circulant_size;
    int n;int m;
    fscanf(fin,"%d",&circulant_size);
    n=circulant_size*bn;m=circulant_size*bm;
    fprintf(fout,"%d %d\n",m,n);
    for(int i=1;i<bn+1;i++){
        for(int j=1;j<bm+1;j++){
            fscanf(fin,"%d",&mask[i][j]);
        }
    }
    int** b_mask=(int**)malloc(sizeof(int*)*(n+1));
    for(int i=0;i<n+1;i++)
        b_mask[i]=(int*)malloc(sizeof(int)*(m+1));
    for(int i=1;i<n+1;i++){
        for(int j=1;j<m+1;j++)
            b_mask[i][j]=0;
    }
    for(int i=1;i<bn+1;i++){
        for(int j=1;j<bm+1;j++){
            int degree=b[i][j]*mask[i][j];
            if(degree>0){
                degree--;
                int **temp=NULL;
                temp=identity(circulant_size);
                right_shift(temp,circulant_size,degree,0);
                for(int k=1;k<circulant_size+1;k++){
                    for(int l=1;l<circulant_size+1;l++){
                        b_mask[(i-1)*circulant_size+k][(j-1)*circulant_size+l]=temp[k][l];
                    }
                }
                for(int k=0;k<circulant_size+1;k++)
                    free(temp[k]);
                free(temp);
            }
        }
    }

    int dn=0;int dm=0;
    for(int i=1;i<n+1;i++){
        dm+=b_mask[i][1];
    }
    for(int j=1;j<m+1;j++){
        dn+=b_mask[1][j];
    }
    fprintf(fout,"%d %d\n",dm,dn);
    for(int i=1;i<m+1;i++)  fprintf(fout,"%d ",dm);
    fprintf(fout,"\n");
    for(int i=1;i<n+1;i++)  fprintf(fout,"%d ",dn);
    fprintf(fout,"\n");
    for(int i=1;i<m+1;i++){
        for(int j=1;j<n+1;j++){
            if(b_mask[j][i]==1)
                fprintf(fout,"%d ",j);
        }
        fprintf(fout,"\n");
    }
    for(int i=1;i<n+1;i++){
        for(int j=1;j<m+1;j++){
            if(b_mask[i][j]==1){
                //printf("%d %d\n",i,j);
                fprintf(fout,"%d ",j);
            }
        }
        fprintf(fout,"\n");
    }

    for(int i=0;i<n+1;i++)
        free(b_mask[i]);
    for(int i=0;i<bn+1;i++){
        free(b[i]);free(mask[i]);
    }
    free(b);free(mask);free(b_mask);
    fclose(fin);fclose(fout);
    return 0;
}
int **identity(int circulant_size){
    int** ident=(int**)malloc(sizeof(int*)*(circulant_size+1));
    for(int i=0;i<circulant_size+1;i++)
        ident[i]=(int*)malloc(sizeof(int)*(circulant_size+1));
    for(int i=1;i<circulant_size+1;i++){
        for(int j=1;j<circulant_size+1;j++){
            if(i==j)    ident[i][j]=1;
            else    ident[i][j]=0;
        }
    }
    return ident;
}
void right_shift(int** matrix,int circulant_size, int degree,int initial_degree){
    for(int i=1;i<circulant_size+1;i++){
        for(int j=1;j<circulant_size+1;j++){
            if((i-1+initial_degree)%circulant_size==j-1)    matrix[i][j]=0;
            if((i-1+initial_degree+degree)%circulant_size==j-1)     matrix[i][j]=1;
        }
    }
    return;
}
