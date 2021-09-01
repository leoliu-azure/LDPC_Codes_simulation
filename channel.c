#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

unsigned long long SEED=287;
unsigned long long RANV;
int RANI=0;

void normal(double *n1,double *n2,double sigma);
double Ranq1();
void map_to_AWGNchannel(double *word,int wordlength);
void map_back(double *word,int wordlength);
void add_noise(double *word, int wordlength, double *n1, double *n2, double sigma);
double var(double L1,double L2);
double chk(double L1,double L2);
double e(){
    return exp(1.0);
}
int sgn(double x);
double delta(double L1,double L2);
double get_Lj(double yj, double Es_N0);
double direct_table(double x);
double linear_table(double x);

void normal(double *n1,double *n2,double sigma){
double s;
double x1,x2;
do{
    x1=Ranq1();
    x2=Ranq1();
    x1=2*x1-1;
    x2=2*x2-1;
    s=x1*x1+x2*x2;
}while(s>=1.0);
*n1=sigma*x1*sqrt(-log(s)*2/s);
*n2=sigma*x2*sqrt(-log(s)*2/s);
return;
}

double Ranq1(){
if(RANI==0){
    RANV=SEED ^ 4101842887655102017LL;
    RANV^=RANV>>21;
    RANV^=RANV<<35;
    RANV^=RANV>>4;
    RANV=RANV*2685821657736338717LL;
    RANI++;
}
RANV^=RANV>>21;
RANV^=RANV<<35;
RANV^=RANV>>4;
return RANV*2685821657736338717LL*5.42101086242752217E-20;
}

void map_to_AWGNchannel(double *word, int wordlength){
int i=1;
while(i<wordlength){
    if(word[i]==0)
        word[i]=1;
    else if(word[i]==1)
        word[i]=-1;
    i++;
}
return;
}
void add_noise(double *word, int wordlength, double *n1, double *n2, double sigma){
int i=1;
while(i<wordlength){
    normal(n1,n2,sigma);
    word[i]+= *n1;
    //word[i+1]+= *n2;//printf("%lf\n",word[i]);
    i+=1;
}
return;
}

void map_back(double *word,int wordlength){
int i=1;
while(i<wordlength){
    if(word[i]>0)
        word[i]=0;
    else
        word[i]=1;
    i++;
}
return;
}

double var(double L1,double L2){
    return L1+L2;
}

double chk(double L1,double L2){
    //return log(cosh((L1+L2)/2))-log(cosh((L1-L2)/2));
    double min;
    if(fabs(L1)>fabs(L2))
        min = fabs(L2);
    else min=fabs(L1);
    //printf("%lf\n",sgn(L1)*sgn(L2)*min);
    return (double)min*sgn(L1)*sgn(L2)+delta(L1,L2);
}

int sgn(double x){
    return (x>0)-(x<0);
}

double delta(double L1,double L2){
    //return 0.0;
    return log(1+pow(e(),-fabs(L1+L2)))-log(1+pow(e(),-fabs(L1-L2)));
    //return direct_table(fabs(L1+L2))-direct_table(fabs(L1-L2));
    //return linear_table(fabs(L1+L2))-linear_table(fabs(L1-L2));
}

double get_Lj(double yj, double Es_N0){
    return (double)4*Es_N0*yj;
}

double direct_table(double x){
    if(x>=0&&x<=0.196)   return 0.65;
    else if(x>0.196&&x<=0.433)  return 0.55;
    else if(x>0.433&&x<=0.71)   return 0.45;
    else if(x>0.71&&x<=1.05)    return 0.35;
    else if(x>1.05&&x<=1.508)   return 0.25;
    else if(x>1.508&&x<=2.252)  return 0.15;
    else if(x>2.252&&x<=4.5)    return 0.05;
    else if(x>4.5) return 0.0;
    return -1;
}

double linear_table(double x){
    if(x>=0&&x<=0.5)   return -x/2+0.7;
    else if(x>0.5&&x<=1.6)  return -x/4+0.575;
    else if(x>1.6&&x<=2.2)   return -x/8+0.375;
    else if(x>2.2&&x<=3.2)    return -x/16+0.2375;
    else if(x>3.2&&x<=4.4)   return -x/32+0.1375;
    else if(x>4.4) return 0.0;
    return -1;
}
