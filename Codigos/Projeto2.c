// Thássio Gabriel Farias dos Santos 140163697 //

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

    time_t t;

struct ARG{

    int argv;
};
    int* randomicoGrama;
    int* randomicoAsfalto;
    double *vetorTreinoGrama[25];
    double *vetorTreinoAsfalto[25];
    double *vetorTesteGrama[25];
    double *vetorTesteAsfalto[25];

double *normaliza(double* vetor){

    float menor = 1000;
    float maior = 0;
    float divisao;
    float conta;
    for(int i = 0; i<536;i++){
        if(menor > vetor[i]){
            menor = vetor[i];
        }
        if(maior < vetor[i]){
            maior = vetor[i];
        }
    }
    divisao = maior - menor;
    for(int i = 0; i<536;i++){
    vetor[i] = (vetor[i] - menor) / divisao;
    }

    return vetor;

}

float *GLCM(int **matriz){
   
    float glcmdireita[256][256];
    float glcmesquerda[256][256];
    float glcmcima[256][256];
    float glcmbaixo[256][256];
    float glcmdireitabaixo[256][256];
    float glcmesquerdabaixo[256][256];
    float glcmesquerdacima[256][256];
    float glcmdireitacima[256][256];

    for (int i = 0;i<256;i++){
        for(int j = 0;j<256;j++){
           glcmdireita[i][j] = 0;
           glcmesquerda[i][j] = 0;
           glcmcima[i][j] = 0;
           glcmbaixo[i][j] = 0;
           glcmdireitabaixo[i][j] = 0;
           glcmdireitacima[i][j] = 0;
           glcmesquerdacima[i][j] = 0;
           glcmesquerdabaixo[i][j] = 0;
        }
    }

    for (int i = 0;i<1025;i++){
        for(int j=0;j<1024;j++){
         glcmdireita[matriz[i][j]][matriz[i][j+1]] = glcmdireita[matriz[i][j]][matriz[i][j+1]] + 1;
        }
    }

    for (int i = 0;i<1025;i++){
        for(int j=1;j<1025;j++){
         glcmesquerda[matriz[i][1025-j]][matriz[i][1024-j]] = glcmesquerda[matriz[i][1025-j]][matriz[i][1024-j]] + 1;
        }
    }

    for (int i = 0;i<1025;i++){
        for(int j=1;j<1025;j++){
         glcmcima[matriz[1025-j][i]][matriz[1024-j][i]] = glcmcima[matriz[1025-j][i]][matriz[1024-j][i]] + 1;
        }
    }
 

    for (int i = 0;i<1025;i++){
        for(int j=0;j<1024;j++){
         glcmbaixo[matriz[j][i]][matriz[j+1][i]] = glcmbaixo[matriz[j][i]][matriz[j+1][i]] + 1;
        }
    }

    for (int i = 0;i<1024;i++){
        for(int j=0;j<1024;j++){
         glcmdireitabaixo[matriz[j][i]][matriz[j+1][i+1]] = glcmdireitabaixo[matriz[j][i]][matriz[j+1][i+1]] + 1;
        }
    }

    for (int i = 1;i<1025;i++){
        for(int j=0;j<1024;j++){
         glcmesquerdabaixo[matriz[j][1025-i]][matriz[j+1][1024-i]] = glcmesquerdabaixo[matriz[j][1025-i]][matriz[j+1][1024-i]] + 1;
        }
    }

    for (int i = 1;i<1025;i++){
        for(int j=1;j<1025;j++){
         glcmesquerdacima[matriz[1025-j][1025-i]][matriz[1024-j][1024-i]] = glcmesquerdacima[matriz[1025-j][1025-i]][matriz[1024-j][1024-i]] + 1;
        }
    }


    for (int i = 0;i<1024;i++){
        for(int j=1;j<1025;j++){
         glcmdireitacima[matriz[1025-j][i]][matriz[1024-j][i+1]] = glcmdireitacima[matriz[1025-j][i]][matriz[1024-j][i+1]] + 1;
        }
    }


    float energia[8];
    float contraste[8];
    float homogeneidade[8];
    float modulo;
    float d,n,m;

    for(int i = 0; i<8;i++){

        energia[i] = 0;
        contraste[i] = 0;
        homogeneidade[i] = 0;
    }

    for(int i = 0; i<256;i++){
        for(int j = 0; j<256; j++){

    if(j > i){
    modulo = -1;
    }
    else{
    modulo = 1;
    }
    d = (i-j);
    n = d*d;
    m = 1 + (d*modulo);

        energia[0] = energia[0] + glcmdireita[i][j]*glcmdireita[i][j];
        energia[1] = energia[1] + glcmesquerda[i][j]*glcmesquerda[i][j];
        energia[2] = energia[2] + glcmcima[i][j]*glcmcima[i][j];
        energia[3] = energia[3] + glcmbaixo[i][j]*glcmbaixo[i][j];
        energia[4] = energia[4] + glcmdireitabaixo[i][j]*glcmdireitabaixo[i][j];
        energia[5] = energia[5] + glcmesquerdabaixo[i][j]*glcmesquerdabaixo[i][j];
        energia[6] = energia[6] + glcmesquerdacima[i][j]*glcmesquerdacima[i][j];
        energia[7] = energia[7] + glcmdireitacima[i][j]*glcmdireitacima[i][j];

        contraste[0] = contraste[0] + n*glcmdireita[i][j];
        contraste[1] = contraste[1] + n*glcmesquerda[i][j];
        contraste[2] = contraste[2] + n*glcmcima[i][j];
        contraste[3] = contraste[3] + n*glcmbaixo[i][j];
        contraste[4] = contraste[4] + n*glcmdireitabaixo[i][j];
        contraste[5] = contraste[5] + n*glcmesquerdabaixo[i][j];
        contraste[6] = contraste[6] + n*glcmesquerdacima[i][j];
        contraste[7] = contraste[7] + n*glcmdireitacima[i][j];

        homogeneidade[0] = homogeneidade[0] + ((float)glcmdireita[i][j]/m);
        homogeneidade[1] = homogeneidade[1] + ((float)glcmesquerda[i][j]/m);
        homogeneidade[2] = homogeneidade[2] + ((float)glcmcima[i][j]/m);
        homogeneidade[3] = homogeneidade[3] + ((float)glcmbaixo[i][j]/m);
        homogeneidade[4] = homogeneidade[4] + ((float)glcmdireitabaixo[i][j]/m);
        homogeneidade[5] = homogeneidade[5] + ((float)glcmesquerdabaixo[i][j]/m);
        homogeneidade[6] = homogeneidade[6] + ((float)glcmesquerdacima[i][j]/m);
        homogeneidade[7] = homogeneidade[7] + ((float)glcmdireitacima[i][j]/m);

        }
    }

    float *retorna = malloc(sizeof(float)*24);

    for(int i=0;i<8;i++){
    retorna[i] = energia[i];
    }
    for(int i=0;i<8;i++){
    retorna[i+8] = contraste[i];
    }
    for(int i=0;i<8;i++){
    retorna[i+16] = homogeneidade[i];
    }

    return retorna;

}

int ILBP(int** imagem,int i,int j){

    int ILBP[9];
    float media;
    float matriz[3][3];
    int  matriztemp[3][3];

    media = (imagem[i-1][j-1] + imagem[i-1][j] + imagem[i-1][j+1] + imagem[i][j-1] +
    imagem[i][j] + imagem[i][j+1] + imagem[i+1][j-1] + imagem[i+1][j] + imagem[i+1][j+1])/9;

    if((imagem[i-1][j-1] - media) >= 0){
    matriztemp[0][0] = 1;
    }else{matriztemp[0][0] = 0;}
    if((imagem[i-1][j] - media) >= 0){
    matriztemp[0][1] = 1;
    }else{matriztemp[0][1] = 0;}
    if((imagem[i-1][j+1] - media) >= 0){
    matriztemp[0][2] = 1;
    }else{matriztemp[0][2] = 0;}
    if((imagem[i][j-1] - media) >= 0){
    matriztemp[1][0] = 1;
    }else{matriztemp[1][0] = 0;}
    if((imagem[i][j] - media) >= 0){
    matriztemp[1][1] = 1;
    }else{matriztemp[1][1] = 0;}
    if((imagem[i][j+1] - media) >= 0){
    matriztemp[1][2] = 1;
    }else{matriztemp[1][2] = 0;}
    if((imagem[i+1][j-1] - media) >= 0){
    matriztemp[2][0] = 1;
    }else{matriztemp[2][0] = 0;}
    if((imagem[i+1][j] - media) >= 0){
    matriztemp[2][1] = 1;
    }else{matriztemp[2][1] = 0;}
    if((imagem[i+1][j+1] - media) >= 0){
    matriztemp[2][2] = 1;
    }else{matriztemp[2][2] = 0;}

    ILBP[0] = matriztemp[0][0];
    ILBP[1] = matriztemp[0][1];
    ILBP[2] = matriztemp[0][2];
    ILBP[3] = matriztemp[1][2];
    ILBP[4] = matriztemp[2][2];
    ILBP[5] = matriztemp[2][1];
    ILBP[6] = matriztemp[2][0];
    ILBP[7] = matriztemp[1][0];
    ILBP[8] = matriztemp[1][1];

    //---------------------------- rotação ---------------------------//

    int menorvalor[9];

    int ILBP1[9];

    ILBP1[0] = ILBP[8];
    ILBP1[1] = ILBP[0];
    ILBP1[2] = ILBP[1];
    ILBP1[3] = ILBP[2];
    ILBP1[4] = ILBP[3];
    ILBP1[5] = ILBP[4];
    ILBP1[6] = ILBP[5];
    ILBP1[7] = ILBP[6];
    ILBP1[8] = ILBP[7];

    int valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[8] = valor;
     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[0] = valor;

    ILBP1[0] = ILBP[7];
    ILBP1[1] = ILBP[8];
    ILBP1[2] = ILBP[0];
    ILBP1[3] = ILBP[1];
    ILBP1[4] = ILBP[2];
    ILBP1[5] = ILBP[3];
    ILBP1[6] = ILBP[4];
    ILBP1[7] = ILBP[5];
    ILBP1[8] = ILBP[6];

     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[1] = valor;

    ILBP1[0] = ILBP[6];
    ILBP1[1] = ILBP[7];
    ILBP1[2] = ILBP[8];
    ILBP1[3] = ILBP[0];
    ILBP1[4] = ILBP[1];
    ILBP1[5] = ILBP[2];
    ILBP1[6] = ILBP[3];
    ILBP1[7] = ILBP[4];
    ILBP1[8] = ILBP[5];

     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[2] = valor;

    ILBP1[0] = ILBP[5];
    ILBP1[1] = ILBP[6];
    ILBP1[2] = ILBP[7];
    ILBP1[3] = ILBP[8];
    ILBP1[4] = ILBP[0];
    ILBP1[5] = ILBP[1];
    ILBP1[6] = ILBP[2];
    ILBP1[7] = ILBP[3];
    ILBP1[8] = ILBP[4];

     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[3] = valor;

    ILBP1[0] = ILBP[4];
    ILBP1[1] = ILBP[5];
    ILBP1[2] = ILBP[6];
    ILBP1[3] = ILBP[7];
    ILBP1[4] = ILBP[8];
    ILBP1[5] = ILBP[0];
    ILBP1[6] = ILBP[1];
    ILBP1[7] = ILBP[2];
    ILBP1[8] = ILBP[3];

     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[4] = valor;

    ILBP1[0] = ILBP[3];
    ILBP1[1] = ILBP[4];
    ILBP1[2] = ILBP[5];
    ILBP1[3] = ILBP[6];
    ILBP1[4] = ILBP[7];
    ILBP1[5] = ILBP[8];
    ILBP1[6] = ILBP[0];
    ILBP1[7] = ILBP[1];
    ILBP1[8] = ILBP[2];

     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[5] = valor;

    ILBP1[0] = ILBP[2];
    ILBP1[1] = ILBP[3];
    ILBP1[2] = ILBP[4];
    ILBP1[3] = ILBP[5];
    ILBP1[4] = ILBP[6];
    ILBP1[5] = ILBP[7];
    ILBP1[6] = ILBP[8];
    ILBP1[7] = ILBP[0];
    ILBP1[8] = ILBP[1];

     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[6] = valor;

    ILBP1[0] = ILBP[1];
    ILBP1[1] = ILBP[2];
    ILBP1[2] = ILBP[3];
    ILBP1[3] = ILBP[4];
    ILBP1[4] = ILBP[5];
    ILBP1[5] = ILBP[6];
    ILBP1[6] = ILBP[7];
    ILBP1[7] = ILBP[8];
    ILBP1[8] = ILBP[0];

     valor = 0;
    for(int i = 0;i < 9;i++){
    if(ILBP1[8-i] == 1){  
    valor = valor + pow(2, i);
    }
    }
    menorvalor[7] = valor;


    int menorILBP = 100000;

    for(int i = 0; i<9; i++){

    if(menorILBP > menorvalor[i]){

        menorILBP = menorvalor[i];
    }

    }

    return menorILBP;

}


double *Matriz(char* nomedoarquivo){

    int histograma[512];
    double *VetorCompleto = malloc(sizeof(double)*536);
    double *VetorNormalizado = malloc(sizeof(double)*536);
    int ILBP_imagem;
    float *GLCM_imagem;
    int i,j;
    FILE *file;

     for(int i=0; i<512; i++){
     histograma[i] = 0;
     }

     printf("%s\n", nomedoarquivo);
     file = fopen(nomedoarquivo, "r");
    int **imagem = (int**)malloc(1025*sizeof(int*));
    for(int i=0; i<1025; i++){
        imagem[i]=(int*)malloc(1025*sizeof(int));
    }
        for(int i=0; i<1025; i++){
          for(int j=0; j<1025; j++){
               fscanf(file, "%d%*c", &imagem[i][j]);         
            }
        }
    fclose(file);

    GLCM_imagem = GLCM(imagem); 

    for(i = 1; i<1024;i++){
    for(j = 1; j<1024;j++){

    ILBP_imagem = ILBP(imagem,i,j);
    histograma[ILBP_imagem] =  histograma[ILBP_imagem] + 1;

    }
    }
    free(imagem);

    for(int i = 0; i<512;i++){

    VetorCompleto[i] = histograma[i];
    }

    for(int i = 0; i<24;i++){

    VetorCompleto[512+i] = GLCM_imagem[i];
    }

    VetorNormalizado = normaliza(VetorCompleto);

    return VetorNormalizado;
    
}

int *GeraNumeroAleatorio(){
    char *nomeasfalto;
    char *stringdovalor;
    stringdovalor = malloc(2*sizeof(char));
    int *inicializa = malloc(50*sizeof(int));

    for (int i = 0; i < 50; i++) {
        inicializa[i] = i+1;
    }
    for (int i = 0; i < 50; i++) {
        int temp = inicializa[i];
        int aleatorio = rand() % 50;
        inicializa[i] = inicializa[aleatorio];
        inicializa[aleatorio] = temp;
    }
    return inicializa;
    }

char *geranomegrama(int aleatorios){

    int i;
    char* nome = malloc(10*sizeof(int));
    char stringnumero[20];
        strcpy(nome, "grass/grass_");
        sprintf(stringnumero, "%d", aleatorios);
        if(aleatorios <10){
        strcat(nome, "0");
        }
        sprintf(stringnumero, "%d", aleatorios);
        strcat(nome, stringnumero);
        strcat(nome, ".txt");
    return nome;

}

char *geranomeasfalto(int aleatorios){

    int i;
    char* nome = malloc(10*sizeof(int));
    char stringnumero[20];
        strcpy(nome, "asphalt/asphalt_");
        sprintf(stringnumero, "%d", aleatorios);
        if(aleatorios <10){
        strcat(nome, "0");
        }
        sprintf(stringnumero, "%d", aleatorios);
        strcat(nome, stringnumero);
        strcat(nome, ".txt");

    return nome;

}


void *rodathreadgr(void *valor){
    int n;
    struct ARG *arg = (struct ARG *)valor;
    n = arg->argv;
    if(n < 25){
    vetorTreinoGrama[n] = Matriz(geranomegrama(randomicoGrama[n]));
    }
    else{
    vetorTesteGrama[n-25] = Matriz(geranomegrama(randomicoGrama[n]));
    }
}

void *rodathreadas(void *valor){
    int n;
    struct ARG *arg = (struct ARG *)valor;
    n = arg->argv;
    if(n < 25){
    vetorTreinoAsfalto[n] = Matriz(geranomeasfalto(randomicoAsfalto[n]));
    }
    else{
    vetorTesteAsfalto[n-25] = Matriz(geranomeasfalto(randomicoAsfalto[n]));
    }
}


int main(){
struct ARG args[50];
pthread_t tid[50];
pthread_attr_t attr[50];
pthread_t tid1[50];
pthread_attr_t attr1[50];
    
    int i;
    float mediaGrama[536];
    float mediaAsfalto[536];

    for(int i=0; i<1; i++){
        vetorTreinoGrama[i]=(double*)malloc(536*sizeof(double));
        vetorTreinoAsfalto[i]=(double*)malloc(536*sizeof(double));
        vetorTesteAsfalto[i]=(double*)malloc(536*sizeof(double));
        vetorTesteGrama[i]=(double*)malloc(536*sizeof(double));
    }
    srand((unsigned) time(&t));


    randomicoGrama = GeraNumeroAleatorio();

    randomicoAsfalto = GeraNumeroAleatorio();

    printf("%d\n", randomicoGrama[0]);
    printf("%d\n", randomicoAsfalto[0]);

for(int i = 0; i<50;i++){
args[i].argv = i;
pthread_attr_init(&attr[i]);
pthread_attr_init(&attr1[i]);
pthread_create(&tid[i],&attr[i],rodathreadgr,&args[i]);
pthread_create(&tid1[i],&attr1[i],rodathreadas,&args[i]);
}

for(int i = 0; i<50;i++){

pthread_join(tid[i], NULL);
}

for(int i = 0; i<50;i++){
pthread_join(tid1[i], NULL);
}

    for(i=0;i<536;i++){
    mediaGrama[i] = 0;
    mediaAsfalto[i] = 0;
    }


    for(int i = 0;i <25;i++){
        for(int j=0;j<536;j++){
            mediaGrama[j] = mediaGrama[j] + vetorTreinoGrama[i][j];
            mediaAsfalto[j] = mediaAsfalto[j] + vetorTreinoAsfalto[i][j];
        }
    }


     for(i=0;i<536;i++){
     mediaGrama[i] = mediaGrama[i]/25;
     mediaAsfalto[i] = mediaAsfalto[i]/25;
     }

     double GramadistanciaGrama[25],GramadistanciaAsfalto[25],AsfaltodistanciaAsfalto[25],AsfaltodistanciaGrama[25];
     double z,x,z1,x1;

    for(int i = 0;i<25;i++){
        GramadistanciaGrama[i] = 0;
        GramadistanciaAsfalto[i] = 0;
        AsfaltodistanciaGrama[i] = 0;
        AsfaltodistanciaAsfalto[i] = 0;
    }

    for(int i = 0;i<25;i++){
         for(int j = 0;j<536;j++){

            z = vetorTesteGrama[i][j] - mediaGrama[j];
            x = vetorTesteGrama[i][j] - mediaAsfalto[j];

            z1 = vetorTesteAsfalto[i][j] - mediaGrama[j];
            x1 = vetorTesteAsfalto[i][j]- mediaAsfalto[j];

            GramadistanciaGrama[i] = GramadistanciaGrama[i] + pow(z,2);
            GramadistanciaAsfalto[i] = GramadistanciaAsfalto[i] + pow(x,2);

            AsfaltodistanciaGrama[i] = AsfaltodistanciaGrama[i] + pow(z1,2);
            AsfaltodistanciaAsfalto[i] = AsfaltodistanciaAsfalto[i] + pow(x1,2);
         }
    }

    double taxadeacerto, acertoGrama = 0,acertoAsfalto = 0,falsaAceitacao= 0,falsaRejeicao=0;

    for(int i = 0;i<25;i++){


     GramadistanciaGrama[i] = sqrt(GramadistanciaGrama[i]);
     GramadistanciaAsfalto[i] = sqrt(GramadistanciaAsfalto[i]);
     AsfaltodistanciaGrama[i] = sqrt(AsfaltodistanciaGrama[i]);
     AsfaltodistanciaAsfalto[i] = sqrt(AsfaltodistanciaAsfalto[i]);

     if(GramadistanciaAsfalto[i] > GramadistanciaGrama[i]){
        acertoGrama++;
        printf("O arquivo Grama %d foi reconhecido como grama\n", randomicoGrama[i+25]);
     }
     else{
        falsaRejeicao++;
        printf("O arquivo Grama %d foi reconhecido como asfalto\n", randomicoGrama[i+25]);

     }
     if(AsfaltodistanciaGrama[i] > AsfaltodistanciaAsfalto[i]){
        acertoAsfalto++;
        printf("O arquivo Asfalto %d foi reconhecido como Asfalto\n", randomicoAsfalto[i+25]);

     }
     else{
        falsaAceitacao++;
        printf("O arquivo Asfalto %d foi reconhecido como Grama\n", randomicoAsfalto[i+25]);
     }

     }

     taxadeacerto = ((acertoGrama + acertoAsfalto) / 50.0)*100.0;
     falsaAceitacao = (falsaAceitacao / 25.0)*100.0;
     falsaRejeicao = (falsaRejeicao / 25.0)*100.0;

     printf("A taxa de acerto e: %lf%%\n", taxadeacerto);
     printf("A taxa de falsa Aceitacao e: %lf%%\n", falsaAceitacao);
     printf("A taxa de falsa Rejeicao e: %lf%%\n", falsaRejeicao);


}   