//
//  main.c
//  anagauss
//
//  Created by Monica  on 08/04/23.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    char filename[] = "/Users/monicacesario/Desktop/dati.txt";
    FILE *fp = fopen(filename, "r");

    if (fp == NULL) {
        printf("Errore: impossibile aprire il file %s.\n", filename);
        exit(EXIT_FAILURE);
    }

    float sum = 0;
    int count = 0;
    float num,avg=0.0,sum2=0.0,var=0.0;
    int num_discarded = 0;

    while (fscanf(fp, "%f", &num) != EOF) {
        if (num_discarded < 5000) {
            num_discarded++;
            continue;
        }
        sum += num;
        sum2 += pow(num-avg,2);
        count++;
    }
    
   
    
   
    
    fclose(fp);
    avg = sum / count;
    var=sum2/(count*(count-1));
    float sigma =sqrt(var);
    printf("Media: %f Deviazione standard: %f\n", avg,sigma);

    return 0;
}
