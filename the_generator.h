#ifndef THE_HEADER_H_INCLUDED
#define THE_HEADER_H_INCLUDED
void generate_information_bits(int wanted_length,int *bit_seed, int *bits_process);
void reduced_row_echelon_form(int **row_echelon,int row_count,int column_count);
void add_or_sub(int **row_echelon,int dest, int source, int scalar, int m);
void row_swap(int **row_echelon,int i,int r,int m);
int bin_add(int a, int b);
int bin_multiply(int a, int b);
int **ano_approach(int **h, int n, int m);

#endif // THE_HEADER_H_INCLUDED
