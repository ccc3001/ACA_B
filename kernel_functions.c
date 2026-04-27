
#include "kernel_functions.h"
#include <math.h>
#include <assert.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 4096
#define MAX(a, b) ((a) > (b) ? (a) : (b))
int count_columns(char* line)
{
    int count = 0;
    char* tmp = strdup(line);
    char* tok = strtok(tmp, ",\n");

    while (tok) {
        count++;
        tok = strtok(NULL, ",\n");
    }

    free(tmp);
    return count;
}

double** read_csv_matrix(const char* filename, int n, int* out_cols)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("Error opening file\n");
        return NULL;
    }

    char line[MAX_LINE];

    // --- Read first line ---
    if (!fgets(line, MAX_LINE, fp)) {
        fclose(fp);
        return NULL;
    }

    int cols = count_columns(line);
    *out_cols = cols;

    // --- Allocate matrix (column-major like you used) ---
    double** matrix = malloc(cols * sizeof(double*));
    for (int j = 0; j < cols; j++)
        matrix[j] = malloc(n * sizeof(double));

    // --- Fill first row ---
    char* tok = strtok(line, ",\n");
    for (int j = 0; j < cols && tok; j++) {
        matrix[j][0] = atof(tok);
        tok = strtok(NULL, ",\n");
    }

    // --- Read remaining rows ---
    int i = 1;
    while (fgets(line, MAX_LINE, fp) && i < n)
    {
        tok = strtok(line, ",\n");
        for (int j = 0; j < cols && tok; j++) {
            matrix[j][i] = atof(tok);
            tok = strtok(NULL, ",\n");
        }
        i++;
    }

    fclose(fp);
    return matrix;
}





pfullmatrix
gaussian_kernel_matrix(int rows, int cols,float h){
    assert(rows >= 0);
    assert(cols >= 0);
    assert(h > 0);

    int i,j,k,SUSY_cols;
    double** matrix;
    matrix=read_csv_matrix("data/susy.csv",MAX(rows,cols),&SUSY_cols);

    double h2=2* h * h;
    pfullmatrix gaussian_matrix = new_fullmatrix(rows,cols); 
    for (i = 0; i< gaussian_matrix->rows;i++){
        for (j=0;j< gaussian_matrix->cols;j++){
            double sum_dist=0;
            for (k=0;k<8;k++){
                double dist = matrix[i][k]-matrix[j][k];
                sum_dist +=dist*dist;
            }
            gaussian_matrix->e[j*gaussian_matrix->rows+i]=exp(sum_dist/h2);
        };
    };
    return gaussian_matrix;
}


int main()
{
    pfullmatrix matrix;
    matrix = gaussian_kernel_matrix(50,50,0.2);
    print_fullmatrix(matrix);
    return 1;
}