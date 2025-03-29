//
//  main.c
//  ising_metro con ran2()
//
//  Created by Monica  on 11/04/23.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


float ran2(void);
void ranstart(void);
void ranfinish(void);

void geometry(int, int*);
void initialize_lattice(int, int, float**);
float energy(float, int, float**, int*);
float magnetization(float, int, float**);
void update_metropolis(int,float,float, int*, float**);
void print_matrix(int, float**);

long idum = 0, idum2 = 0, iy = 0, iv[NTAB];

int main(void) {
    const int nlatt=40;
    float xmagn=0.0,xene=0.0;
    
    double xmagn_sum = 0.0, xene_sum = 0.0;

    int iflag = 0, measures , i_decorrel ;
    double extfield,beta;
    
    float **field = (float **)malloc(nlatt * sizeof(float*));
    for(int i = 0; i < nlatt; i++) field[i] = (float *)malloc(nlatt * sizeof(float));

    int val[2*nlatt];
    


    FILE *input_file, *lattice_file, *data_file;

    // apertura del file da cui leggere i parametri della simulazione
    input_file = fopen("/Users/monicacesario/Desktop/modulo1/input.txt", "r");

    if (input_file == NULL) {
        printf("Errore: il file input non può essere aperto.\n");
        return 1;
    }

    printf("PARAMETERS\n");

    fscanf(input_file, "%d", &iflag);
    printf("IFLAG: %d\n", iflag);
    fscanf(input_file, "%d", &measures);
    printf("MEASURES: %d\n", measures);
    fscanf(input_file, "%d", &i_decorrel);
    printf("DECORREL: %d\n", i_decorrel);
    fscanf(input_file, "%lf", &extfield);
    printf("EXTFIELD: %lf\n", extfield);
    fscanf(input_file, "%lf", &beta);
    printf("BETA: %lf\n", beta);

    printf("\n\n");

    // apertura di un file sul quale scrivere le misure della magnetizzazione
    
    data_file = fopen("/Users/monicacesario/Desktop/modulo1/isto40_0.44.dat", "w");
    if (data_file == NULL) {
        printf("Errore nell'apertura del file data!\n");
        return 1;
    }

    // inizializzo generatore di numeri random
    ranstart();

   // for (beta = 0.441; beta < 0.83; beta += 0.0125){
        
        // inizializzo condizioni al bordo
        geometry(nlatt, (int *)val);
        
        // inizializzo configurazione iniziale
        initialize_lattice(iflag, nlatt, (float **)field);
        
       
       
            // misura delle osservabili fisiche
            xmagn = magnetization(xmagn, nlatt, (float **)field);
            xene = energy(xene, nlatt, (float **)field, (int *)val);
            printf("%f %f\n",xmagn, xene);
            
            for (int i = 0; i < measures; i++) {
                
                update_metropolis(nlatt, beta, extfield, (int *)val, (float **)field);
                
                if ((i + 1) % i_decorrel == 0) {
                    // misura delle osservabili fisiche
                    xmagn = (magnetization(xmagn, nlatt, (float **)field));
                    xene = energy(xene, nlatt, (float **)field, (int *)val);
                    xmagn_sum += xmagn;
                    xene_sum += xene;
                    printf("%d %f %f %f\n", i + 1, beta, xmagn, xene);
                    // scrivo misure sul file data per poi fare l'analisi
                    fprintf(data_file, "%d %f %f %f\n", i,beta, xmagn, xene);
                    fflush(data_file);
                }
         //   }
            
            fprintf(data_file, "\n");
        }
        
  

    // salvo la configurazione per poter eventualmente ripartire
    lattice_file = fopen("/Users/monicacesario/Desktop/modulo1/lattice_25.txt", "w");
    if (lattice_file == NULL) {
        printf("Errore: il file lattice non può essere aperto.\n");
        return 1;
    }
    for (int i = 0; i < nlatt; i++) {
        for (int j = 0; j < nlatt; j++) {
            fprintf(lattice_file, "%f ", field[i][j]);
        }
        fprintf(lattice_file, "\n");
    }

    ranfinish();

    fclose(input_file);
    fclose(data_file);
    fclose(lattice_file);

    printf("\n FINAL \n");
    print_matrix(nlatt, (float **)field);

    printf("\nDONE!\n");
    return 0;
    }


// Per ogni coordinata definisco il passo in avanti o indietro con le opportune condizioni al bordo
void geometry(int nlatt, int* val) {
  
    for(int i=0; i<nlatt-1; ++i){
        val[i]=i+1;
    }
    val[nlatt-1]=0;

    for(int i=nlatt+1; i<2*nlatt; ++i){
        val[i]=(i-1)-nlatt;
    }
    val[nlatt]=nlatt-1;
       
}
//assegno la configurazione di partenza alla catena di Markow
void initialize_lattice(int iflag, int nlatt, float **field) {
  
    double x;
    
    FILE *lattice_file;
   
    // matrice con le variabili di spin
    
    lattice_file = fopen("/Users/monicacesario/Desktop/modulo1/lattice40_0.44.txt", "r");
    if (lattice_file == NULL) {
        printf("Errore: il file lattice non può essere aperto.\n");
        
    }
   
    //partenza a freddo (tutti gli spin a 1 come se fosse T = 0)
    if(iflag <= 0){
        for(int i = 0; i<nlatt; i++){
            for(int j = 0; j<nlatt; j++){
                field[i][j]=1.0;
            }
        }
    }
    
    //partenza a caldo ( spin random come e fosse T = infinito)
    else if(iflag == 1){
        for(int i=0; i<nlatt; i++){
            for(int j=0; j<nlatt; j++){
                //nro random tra 0 e 1
                x = (double)ran2();
                field[i][j]=1.0;
                if(x <= 0.5) field[i][j]=-1.0;
            }
        }
        
    }
    
   // da dove ero rimasto l'ultima volta, devo leggere dal file lattice i vlori di field e seed
    else {
        for (int i=0; i<nlatt; i++){
            for(int j=0; j<nlatt; j++){
                fscanf(lattice_file, "%f ", &field[i][j]);
            }
        }
       
        }
       
     
        fclose(lattice_file);

    }
    
    //calcolo la magnetizzazione del reticolo
    
    float magnetization(float xmagn, int nlatt, float **field){
        
        int nvol = pow(nlatt,2);
        // inizializzo magnetizzazione
        xmagn=0.0;
        
        //loop su tutto il reticolo
        for(int i = 0; i<nlatt; i++){
            for(int j= 0; j<nlatt; j++){
                //e sommo tutti i valri del campo
                xmagn += field[i][j];
            }
        }
        // normalizzo dividendo per il volume
        return  xmagn/(float)nvol;
        
        
    }
    
    
    // energia media (=0 per configurazione ordinata e campo esterno 0)
    float energy(float xene, int nlatt, float **field, int *val){
        
        int nvol = pow(nlatt,2);
        int npp[nlatt],nmm[nlatt];
        int ip,im,jp,jm;
        float force=0.0,extfield=0.0;
        
       
         for(int i = 0; i < nlatt; i++){
            npp[i]=val[i];
        }

        for(int i=nlatt; i< 2*nlatt; i++){
            nmm[i-nlatt] = val[i];
        }
      
       // loop sul reticolo
        for(int i = 0; i<nlatt; i++){
            for(int j=0; j<nlatt; j++){
                //calcolo coordinate prime vicine del sito
                ip=npp[i];
                im=nmm[i];
                jp=npp[j];
                jm=nmm[j];
                
                //somma dei 4 primi vicini
                force = field[i][jp]+field[i][jm]+
                field[ip][j]+field[im][j];
                
                xene = xene - 0.5*force*field[i][j];
                
                //contributo campo esterno
                xene = xene - extfield*field[i][j];
            }
        }
       
        return xene/(float)nvol;
        
    }
    

    void update_metropolis(int nlatt,float beta,float extfield, int *val, float **field){
       
        int npp[nlatt], nmm[nlatt];
        int ip,im,jp,jm,l,m;
        float force,phi;
        double p_rat;
        double x,y,z;
        
        
        for(int i = 0; i < nlatt; i++){
            npp[i]=val[i];
        }

        for(int i=nlatt; i < 2*nlatt; i++){
            nmm[i-nlatt] = val[i];
        }
       
        //loop su tutti i siti
        for (int i = 0; i<nlatt; i++){
            for ( int j= 0; j<nlatt; j++){
                //genero 2 nri random tra 0 e 1
               
                x = (double)ran2();
                y = (double)ran2();
                
                if(x==1.0){
                    printf("x è uno\n");
                } else if (y== 1.0){
                    printf("y è uno\n");
                }
                
                
                l = (int)(x*nlatt);
                m = (int)(y*nlatt);
               
                //calcolo le coordinate dei 4 primi vicini del sito che ho scelto
                
                ip = npp[l];
                im = nmm[l];
                jp = npp[m];
                jm = nmm[m];
                // calcolo la somma dei primi quattro vicini
                
                force = field[l][jp]+field[l][jm]+field[ip][m]+field[im][m];
               
                force = beta*(force + extfield);
                
                // valore attuale dello spin
                phi = field[l][m];
                
                // il tentativo è sempre invertire lo spin. calcolo il rapporto p_rat di prob tra il caso invertito e non
                p_rat = exp(-2.0*phi*force);
               
                // METRO-TEST
                z=(double)ran2();
                
                // x<p_rat verifica anche il caso p_rat > 1, se sì accetto
                if(z<= p_rat){
                    field[l][m]=-phi;
                    
                }
                
                
            }
        }
       
        
    }
    

    
void print_matrix(int nlatt, float* field[]){

    for(int i=0; i<nlatt; i++) {
    // inner loop for column
        for(int j=0; j<nlatt; j++) {
            printf("%lf ", field[i][j]);
        }
        printf("\n"); // new line
    }
}


float ran2(void)

{

    int j;

    long k;

    float temp;

   
    if (idum <= 0) {

        if (-(idum) < 1) idum=1;

        else idum = -(idum);

        idum2=(idum);

        for (j=NTAB+7;j>=0;j--) {

            k=(idum)/IQ1;

            idum=IA1*(idum-k*IQ1)-k*IR1;

            if (idum < 0) idum += IM1;

            if (j < NTAB) iv[j] = idum;

        }

        iy=iv[0];

    }

    k=(idum)/IQ1;

    idum=IA1*(idum-k*IQ1)-k*IR1;

    if (idum < 0) idum += IM1;

    k=idum2/IQ2;

    idum2=IA2*(idum2-k*IQ2)-k*IR2;

    if (idum2 < 0) idum2 += IM2;

    j=iy/NDIV;

    iy=iv[j]-idum2;

    iv[j] = idum;

    if (iy < 1) iy += IMM1;

    if ((temp=AM*iy) > RNMX) return RNMX;

    else return temp;

}

void ranstart(void) {
    FILE *fp;
    
    int i;
    
    fp = fopen("/Users/monicacesario/Desktop/modulo1/randomseed.txt", "r");
    if (fp == NULL) {
        fprintf(stderr, "Errore nell'apertura del file randomseed\n");
        exit(1);
    }
    
    fscanf(fp, "%li", &idum);
    //printf("idum=%li\n",idum);
    fscanf(fp, "%li", &idum2);
    //printf("idum2=%li\n",idum2);
    for (i = 0; i < 32; i++) {
        fscanf(fp, "%li", &iv[i]);
        //printf("iv=%li\n",iv[i]);
    }
    fscanf(fp, "%li", &iy);
    //printf("iy=%li\n",iy);
    
    if (idum >= 0) {
        idum = -idum - 1;
    }
    
    fclose(fp);
}


void ranfinish(void) {
    FILE *fp;
    
    
    fp = fopen("/Users/monicacesario/Desktop/modulo1/randomseed.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "Errore nell'apertura del file randomseed\n");
        exit(1);
    }
    
    fprintf(fp, "%li\n", idum);
    fprintf(fp, "%li\n", idum2);
    for (int i = 0; i < 32; i++) {
        fprintf(fp, "%li\n", iv[i]);
    }
    fprintf(fp, "%li\n", iy);
    
    fclose(fp);
}
