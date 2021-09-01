#ifndef CHANNEL_H_INCLUDED
#define CHANNEL_H_INCLUDED
void normal(double *n1,double *n2,double sigma);
double Ranq1();
void map_to_AWGNchannel(double *word,int wordlength);
void map_back(double *word,int wordlength);
void add_noise(double *word, int wordlength, double *n1, double *n2, double sigma);
double var(double L1,double L2);
double chk(double L1,double L2);
double e();
int sgn(double x);
double delta(double L1,double L2);
double get_Lj(double yj, double Es_N0);
double direct_table(double x);
double linear_table(double x);


#endif // CHANNEL_H_INCLUDED
