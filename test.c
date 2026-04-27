#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 4096

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

int main()
{
    FILE *fp = fopen("data/susy.csv", "r");
    if (!fp) {
        printf("Error opening file\n");
        return 1;
    }

    int n;
    printf("Number of rows: ");
    scanf("%d", &n);

    char line[MAX_LINE];

    // --- Step 1: read first line to determine columns ---
    fgets(line, MAX_LINE, fp);
    int cols = count_columns(line);

    // --- Step 2: allocate matrix ---
    double** matrix = malloc(cols * sizeof(double*));
    for (int i = 0; i < cols; i++)
        matrix[i] = malloc(n * sizeof(double));

    // --- Step 3: process first line ---
    char* tok = strtok(line, ",\n");
    for (int j = 0; j < cols && tok; j++) {
        matrix[j][0] = atof(tok);
        tok = strtok(NULL, ",\n");
    }

    // --- Step 4: read remaining rows ---
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

    // --- Example: print matrix ---
    for (int r = 0; r < cols; r++) {
        for (int c = 0; c < n; c++) {
            printf("%f ", matrix[c][r]);
        }
        printf("\n");
    }

    // --- Free memory ---
    for (int i = 0; i < cols; i++)
        free(matrix[i]);
    free(matrix);

    return 0;
}