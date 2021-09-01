#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "the_generator.h"
#include "channel.h"
#include "vector_decode.h"

int main(const int argc,const char* argv[])
{
    FILE *fin=fopen(argv[1],"r");
    FILE *fout=fopen(argv[2],"w+");

    if(fin==NULL){
        perror("error when opening the file.\n");
        return(-1);
    }

    int m,n,dm,dn;
    int check1=0,check2=0;
    fscanf(fin,"%d %d %d %d",&m,&n,&dm,&dn);
    for(int i=0;i<m;i++){
        int temp;
        fscanf(fin,"%d",&temp);
        check1+=temp;
    }
    for(int i=0;i<n;i++){
        int temp;
        fscanf(fin,"%d",&temp);
        check2+=temp;
    }
    if(check1==check2)  ;
    else{
        perror("wrong input!\n");
        return(-1);
    }

    int **h=(int**)malloc(sizeof(int*)*(n+1));
    int **h_prime=(int**)malloc(sizeof(int*)*(n+1));
    int **row_ref=(int**)malloc(sizeof(int*)*(n+1));
    int **col_ref=(int**)malloc(sizeof(int*)*(m+1));

    for(int i=0;i<n+1;i++){
            h[i]=(int*)malloc(sizeof(int)*(m+1));
            h_prime[i]=(int*)malloc(sizeof(int)*(m+1));
            row_ref[i]=(int*)malloc(sizeof(int)*(dn+1));
    }
    for(int i=0;i<m+1;i++)  col_ref[i]=(int*)malloc(sizeof(int)*(dm+1));

    for(int i=0;i<n+1;i++){
        for(int j=0;j<m+1;j++){
            h[i][j]=0;
            h_prime[i][j]=0;
        }
    }
    for(int j=1;j<m+1;j++){
        for(int i=1;i<dm+1;i++){
            int temp;
            fscanf(fin,"%d",&temp);
            col_ref[j][i]=temp;
            h[temp][j]=1;
            h_prime[temp][j]=1;
        }
    }
    for(int i=1;i<n+1;i++){
        for(int j=1;j<dn+1;j++) {
            int temp;
            fscanf(fin,"%d",&temp);
            row_ref[i][j]=temp;
        }
    }
    reduced_row_echelon_form(h_prime,n,m);printf("here1\n");
    int **g=ano_approach(h_prime,n,m);printf("here2\n");

    double *n1=malloc(sizeof(double)),*n2=malloc(sizeof(double));
    *n1=0,*n2=0;
    double Eb_N0;
    printf("Please input Eb over N0\n");
    scanf("%lf",&Eb_N0);
    double R,sigma;
    R=(double)dm/dn;//doubtful!!!
    Eb_N0=pow(10,(double)(Eb_N0/10));
    sigma=sqrt((double)0.5/R/Eb_N0);

    int max_wrong;
    printf("Input max wrong number:\n");
    scanf("%d",&max_wrong);

    int wrong=0;
    int total=0;
    int cut=0;

    int* bit_seed_init=(int*)malloc(sizeof(int)*6);
    memset(bit_seed_init,0,sizeof(int)*6);
    bit_seed_init[0]=1;

    while(wrong<max_wrong){
        int pass=0;
        double *word=(double*)malloc(sizeof(double)*(m+1));
        memset(word,0,sizeof(double)*(m+1));
        int *information_bit=(int *)malloc(sizeof(int)*(m-n+1));
        generate_information_bits(m-n,bit_seed_init,information_bit);
        for(int i=1;i<m+1;i++){
            for(int j=1;j<m-n+1;j++){
                word[i]=bin_add(word[i],bin_multiply(g[j][i],information_bit[j]));
            }
        }
        int *word_copy=(int*)malloc(sizeof(int)*(m+1));
        for(int i=1;i<m+1;i++)  word_copy[i]=(int)word[i];

        map_to_AWGNchannel(word,m+1);

        add_noise(word,m+1,n1,n2,sigma);

        double *Lj=(double*)malloc(sizeof(double)*(m+1));
        for(int j=1;j<m+1;j++){
            Lj[j]=get_Lj(word[j],R*Eb_N0);
        }

        //Initialization
        double **q_arr=(double**)malloc(sizeof(double)*(n+1));
        for (int i=0;i<n+1;i++){
            q_arr[i]=(double*)malloc(sizeof(double)*(m+1));
        }
        double **u_arr=(double**)malloc(sizeof(double)*(n+1));
        for (int i=0;i<n+1;i++){
            u_arr[i]=(double*)malloc(sizeof(double)*(m+1));
        }
        double **q_store=(double**)malloc(sizeof(double)*(n+1));
        for (int i=0;i<n+1;i++){
            q_store[i]=(double*)malloc(sizeof(double)*(m+1));
        }
        double **u_store=(double**)malloc(sizeof(double)*(n+1));
        for (int i=0;i<n+1;i++){
            u_store[i]=(double*)malloc(sizeof(double)*(m+1));
        }

        for(int i=0;i<n+1;i++){
            for(int j=0;j<m+1;j++){
                q_arr[i][j]=0;
                u_arr[i][j]=0;
                q_store[i][j]=0;
                u_store[i][j]=0;
            }
        }
        for(int i=1;i<n+1;i++){
            for(int j=1;j<m+1;j++){
                if(h[i][j]==1){
                    q_store[i][j]=Lj[j];
                }
            }
        }
        int iteration=0;int error_num=0;int preset_max_iter=80;int flag=0;
        while(flag==0){
            if(iteration>preset_max_iter) break;
            error_num=0;
        //message-passing
        //1.horizontal
            for(int i=1;i<n+1;i++){
                /*int *count=(int*)malloc(sizeof(int)*(m+1));memset(count,0,(m+1)*sizeof(int));
                int num=0;
                for(int j=1;j<m+1;j++){
                    if(h[i][j]==1){
                        count[num]=j;
                        num++;
                    }
                }
                for(int j=1;j<m+1;j++){
                    if(h[i][j]==1){
                            int u=0;
                        for(int k=0;k<dn;k++){
                            if(count[k]==j)
                                continue;
                            else if(u==0){
                                u_arr[i][j]=q_store[i][count[k]];
                                u++;
                            }
                            else{
                                u_arr[i][j]=chk(u_arr[i][j],q_store[i][count[k]]);
                                u++;
                            }
                        }
                    }
                }
                free(count);*/
                for(int j=1;j<dn+1;j++){
                    int u=0;
                    for(int k=1;k<dn+1;k++){
                        if(j==k)    continue;
                        if(u==0){
                            u_arr[i][row_ref[i][j]]=q_store[i][row_ref[i][k]];
                            u++;
                        }
                        else{
                            u_arr[i][row_ref[i][j]]=chk(u_arr[i][row_ref[i][j]],q_store[i][row_ref[i][k]]);
                            u++;
                        }
                    }
                }
            }
            for(int i=1;i<n+1;i++){
                for(int j=1;j<dn+1;j++){
                    u_store[i][row_ref[i][j]]=u_arr[i][row_ref[i][j]];
                    u_arr[i][row_ref[i][j]]=0;
                }
            }

        //2.vertical
            for(int j=1;j<m+1;j++){
                /*int *count=(int*)malloc(sizeof(int)*(n+1));memset(count,0,(n+1)*sizeof(int));
                int num=0;
                for(int i=1;i<n+1;i++){
                    if(h[i][j]==1){
                        count[num]=i;
                        num++;
                    }
                }
                for(int i=1;i<n+1;i++){
                    if(h[i][j]==1){
                        for(int k=0;k<dm;k++){
                            if(count[k]==i)
                                continue;
                            else{
                                q_arr[i][j]=var(q_arr[i][j],u_store[count[k]][j]);
                            }
                        }
                        q_arr[i][j]=var(q_arr[i][j],Lj[j]);
                    }
                }
                free(count);*/
                for(int i=1;i<dm+1;i++){
                    for(int k=1;k<dm+1;k++){
                        if(i==k)    continue;
                        q_arr[col_ref[j][i]][j]=var(q_arr[col_ref[j][i]][j],u_store[col_ref[j][k]][j]);
                    }
                    q_arr[col_ref[j][i]][j]=var(q_arr[col_ref[j][i]][j],Lj[j]);
                }
            }
            for(int i=1;i<dm+1;i++){
                for(int j=1;j<m+1;j++){
                    q_store[col_ref[j][i]][j]=q_arr[col_ref[j][i]][j];
                    q_arr[col_ref[j][i]][j]=0;
                }
            }

        //3.decision
            double *qj=(double*)malloc(sizeof(double)*(m+1));memset(qj,0,sizeof(double)*(m+1));
            int *x_bar=(int*)malloc(sizeof(int)*(m+1));memset(x_bar,0,sizeof(int)*(m+1));

            for(int j=1;j<m+1;j++){
                for(int i=1;i<dm+1;i++){
                        qj[j]=var(qj[j],u_store[col_ref[j][i]][j]);
                }
                qj[j]=var(qj[j],Lj[j]);
            }

            for(int j=1;j<m+1;j++){
                if(qj[j]<0)
                    x_bar[j]=1;
                else
                    x_bar[j]=0;
            }
            for(int i=1;i<m+1;i++){
                if(word_copy[i]!=x_bar[i]){
                    error_num++;
                }
            }
            int parity_check=0;
            int* row_parity_sum=malloc(sizeof(int)*(m+1));
            memset(row_parity_sum,0,sizeof(int)*(m+1));
            for(int i=1;i<n+1;i++){
                for(int j=1;j<m+1;j++){
                    row_parity_sum[i]=bin_add(row_parity_sum[i],bin_multiply(h[i][j],x_bar[j]));
                }
                parity_check+=row_parity_sum[i];
            }
            if(parity_check==0)
                    flag=1;
            pass=error_num;
            free(qj);free(x_bar);free(row_parity_sum);
            iteration++;
        }
    wrong+=pass;
    total+=1;
    for(int i=0;i<n+1;i++){
        free(q_arr[i]);free(u_arr[i]);free(u_store[i]);free(q_store[i]);
    }
        free(q_arr);free(u_arr);free(u_store);free(q_store);
        free(word);free(word_copy);free(information_bit);
        free(Lj);
        if(wrong>cut){
         printf("%d %d\n",total,wrong);
         cut+=500;
        }
    }

    double bit_err = (double)wrong/total/m;
    printf("final result %lf\n",bit_err);

    for(int i=0;i<n+1;i++){
        free(h[i]);
        free(h_prime[i]);
        free(row_ref[i]);
    }
    for(int i=0;i<m+1;i++)  free(col_ref[i]);
    for(int i=0;i<m-n+1;i++)    free(g[i]);
    free(h);free(h_prime);free(row_ref);free(col_ref);
    free(g);
    free(n1);free(n2);
    fclose(fin);fclose(fout);
    return 0;
}

