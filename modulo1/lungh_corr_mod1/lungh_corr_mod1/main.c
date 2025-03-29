/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double calculate_PPx(int r, double ave, double** grid, int dimx, int dimy) {
    double result = 0.0;
    
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            double e1 = grid[i][j];
            double e2 = grid[(i + r) % dimx][j];
            result += (e1 * e2) - ave;
        }
    }
    
    return result / (dimx * dimy);
}

double calculate_PPy(int r, double ave, double** grid, int dimx, int dimy) {
    double result = 0.0;
    
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            double e1 = grid[i][j];
            double e2 = grid[i][(j + r) % dimy];
            result += (e1 * e2) - ave;
        }
    }
    
    return result / (dimx * dimy);
}

int main(void) {
    int r = 1;
    double ave = 0.0;
    int dimx = 10;
    int dimy = 10;
    
    // Esempio di creazione e inizializzazione della griglia
    double** grid = (double**)malloc(dimx * sizeof(double*));
    for (int i = 0; i < dimx; i++) {
        grid[i] = (double*)malloc(dimy * sizeof(double));
        for (int j = 0; j < dimy; j++) {
            grid[i][j] = i * j; // Esempio di inizializzazione della griglia
        }
    }
    
    // Calcolo di PPx
    double PPx = calculate_PPx(r, ave, grid, dimx, dimy);
    printf("PPx: %lf\n", PPx);
    
    // Calcolo di PPy
    double PPy = calculate_PPy(r, ave, grid, dimx, dimy);
    printf("PPy: %lf\n", PPy);
    
    // Deallocazione della griglia
    for (int i = 0; i < dimx; i++) {
        free(grid[i]);
    }
    free(grid);
    
    return 0;
}

*/

#include <stdio.h>

double calculate_correlation(int* grid, int dimx, int dimy, int r) {
    double e1 = 0.0;
    double e2 = 0.0;
    double corr = 0.0;

    // Calcola il valore di e1 come la media degli elementi della griglia
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            e1 += grid[i * dimy + j];
        }
    }
    e1 /= (dimx * dimy);

    // Calcola il valore di e2 come la media degli elementi della griglia "ruotata" di r posizioni
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            int x = (i + r) % dimx;
            e2 += grid[x * dimy + j];
        }
    }
    e2 /= (dimx * dimy);

    // Calcola la correlazione
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            corr += (grid[i * dimy + j] - e1) * (grid[i * dimy + j] - e2);
        }
    }
    corr /= (dimx * dimy);

    return corr;
}

int main() {
    // Esempio di utilizzo della funzione calculate_correlation

    int grid[4][4] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12},
        {13, 14, 15, 16}
    };

    int dimx = 4;
    int dimy = 4;
    int r = 1;

    double corr = calculate_correlation((int*)grid, dimx, dimy, r);
    printf("Correlation: %lf\n", corr);

    return 0;
}

