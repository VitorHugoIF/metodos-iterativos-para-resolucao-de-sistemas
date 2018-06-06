#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int maxN=100; //Limite para o numero de variaveis

int lerEliminacaoGauss(double mat[maxN][maxN+1])
{
    ifstream arq("eliminacaoGauss.txt"); //abre o arquivo matriz.txt para leitura
    if (!arq.is_open())         //teste de abertura
    {
        cout<<"Erro na abertura do arquivo."<<endl;
        exit(1);
    }
    int N;
    arq>>N;
    for( int i=0; i<N; i++)
        for( int j=0; j<N+1; j++)
            arq>>mat[i][j];

    for( int i=0; i<N; i++)
    {
        for( int j=0; j<N+1; j++)
        {
            printf("%7g| ", mat[i][j]);
        }
        cout<<endl;
    }
    cout<<"_______________________________________________"<<endl;
    arq.close();
    return N;

}
int lerJacobi(double mat[maxN][maxN+1])
{
    ifstream arq("jacobi.txt"); //abre o arquivo matriz.txt para leitura
    if (!arq.is_open())         //teste de abertura
    {
        cout<<"Erro na abertura do arquivo."<<endl;
        exit(1);
    }
    int N;
    arq>>N;
    for( int i=0; i<N; i++)
        for( int j=0; j<N+1; j++)
            arq>>mat[i][j];

    for( int i=0; i<N; i++)
    {
        for( int j=0; j<N+1; j++)
        {
            printf("%7g| ", mat[i][j]);
        }
        cout<<endl;
    }
    cout<<"_______________________________________________"<<endl;
    arq.close();
    return N;

}
int lerSeidel(double mat[maxN][maxN+1])
{
    ifstream arq("seidel.txt"); //abre o arquivo matriz.txt para leitura
    if (!arq.is_open())         //teste de abertura
    {
        cout<<"Erro na abertura do arquivo."<<endl;
        exit(1);
    }
    int N;
    arq>>N;
    for( int i=0; i<N; i++)
        for( int j=0; j<N+1; j++)
            arq>>mat[i][j];

    for( int i=0; i<N; i++)
    {
        for( int j=0; j<N+1; j++)
        {
            printf("%7g| ", mat[i][j]);
        }
        cout<<endl;
    }
    cout<<"_______________________________________________"<<endl;
    arq.close();
    return N;

}
bool linhas(double mat[maxN][maxN+1], int tam)
{
    double maximo = 0;
    double soma = 0;

    for( int i=0; i<tam; i++)
    {
        for( int j=0; j<tam; j++)
        {
            if(i!=j)
            {
                soma+= fabs(mat[i][j]);
            }

        }
        soma = soma/mat[i][i];
        if(soma > maximo)
        {
            maximo = soma;
        }
    }
    if(maximo < 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}
bool colunas(double mat[maxN][maxN+1], int tam)
{
    double maximo = 0;
    double soma = 0;

    for( int j=0; j<tam; j++)
    {
        for( int i=0; i<tam; i++)
        {
            if(i!=j)
            {
                soma+= fabs(mat[i][j]);
            }

        }
        soma = soma/mat[j][j];
        if(soma > maximo)
        {
            maximo = soma;
        }
    }
    if(maximo < 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}
bool sassenfeld(double mat[maxN][maxN+1], int tam)
{
    double vetBeta[tam];
    for(int i =0; i< tam; i++)
        vetBeta[i]=0;

    for(int j =1; j<tam; j++)
    {
        vetBeta[0]+=fabs(mat[0][j])/fabs(mat[0][0]);
    }

    for(int i =1; i< tam; i++)
    {
        for(int j =0; j<tam; j++)
        {
            while(j<=i-1)
            {
                vetBeta[i] += fabs((mat[i][j])*vetBeta[j])/fabs(mat[i][i]);
                j++;
            }
            j=i+1;
            vetBeta[i] = vetBeta[i]+ (fabs(mat[i][j])/fabs(mat[i][i]));
        }
    }
    double maior = vetBeta[0];
    for(int i =1; i<tam; i++)
        if(vetBeta[i]>maior)
            maior=vetBeta[i];
    if(maior < 1)
        return true;
    else
        return false;


}
void condiciona(double mat[maxN][maxN+1], int tam)
{
    for( int i=0; i<tam; i++)
    {
        double aii = mat[i][i];
        for( int j=0; j<tam+1; j++)
        {
            mat[i][j] = mat[i][j]/aii;
        }
    }
    for( int i=0; i<tam; i++)
    {
        for( int j=0; j<tam+1; j++)
        {
            printf("%7g| ", mat[i][j]);
        }
        cout<<endl;
    }
    cout<<"_______________________________________________"<<endl;

}
double erro(double *ant, double *atual, int tam)
{
    double aux[tam];
    double numerador,denominador;
    for(int i=0; i<tam; i++)
    {
        aux[i] = fabs(atual[i]-ant[i]);
    }
    numerador = fabs(aux[0]);
    denominador=fabs(atual[0]);
    for(int j =0; j<tam; j++)
    {
        if(fabs(atual[j])>denominador)
        {
            denominador = fabs(atual[j]);
        }
        if(aux[j]>numerador)
        {
            numerador = fabs(aux[j]);
        }

    }
    return numerador/denominador;
}
void jacobi(double *inicial,double mat[maxN][maxN+1], int tam, int iteracoes, double epsilon)
{
    int cont=1,k=0;
    double soma=0;
    double solucao_1[tam];
    for(int i; i<tam; i++)
        solucao_1[i]=0;
    cout<<"{";
    for(int i; i<tam; i++)
    {
        cout<<""<<inicial[i];
        cout<<",";
    }
    cout<<"}";
    cout<<endl;
    cout<<"_______________________________________________"<<endl;

    if(linhas(mat,tam)||colunas(mat,tam))
    {
        condiciona(mat,tam);
        for(int i=0; i<tam; i++)
        {
            soma = mat[i][tam];
            for(int j=tam-1; j >= 0; j--)
            {
                if(j!=i)
                {
                    soma = soma - (mat[i][j])*inicial[j];
                }
            }
            solucao_1[i] = soma/(mat[i][i]);
            printf("Solucao(%d): X%d = %.5f   ",k,i+1,solucao_1[i]);
            if(cont == 3)
            {
                double err = erro(inicial, solucao_1,tam);
                printf("Erro: %.5f  ",err);
                cout<<endl;
                cont=0;
            }
            cont++;
        }
        k++;

        do
        {
            for(int i=0; i<tam; i++)
                inicial[i]=solucao_1[i];

            for(int i=0; i<tam; i++)
            {
                soma = mat[i][tam];
                for(int j=tam-1; j >= 0; j--)
                {
                    if(j!=i)
                    {
                        soma = soma - (mat[i][j])*inicial[j];
                    }
                }

                solucao_1[i] = soma/(mat[i][i]);
                printf("Solucao(%d): X%d = %.5f   ",k,i+1,solucao_1[i]);

                if(cont == 3)
                {
                    double err = erro(inicial, solucao_1,tam);
                    printf("Erro: %.5f  ",err);
                    cout<<endl;
                    cont=0;
                }
                cont++;
            }
            k++;

        }
        while( erro(inicial, solucao_1,tam)>epsilon && k<iteracoes);
        cout<<endl;
        printf("\nSolucoes:\n");
        for(int i=0; i<tam; i++)
            printf("Solucao %d = %.5f ",i+1,solucao_1[i]);
        cout<<endl;
    }
    else
    {
        cout<<"Não atende aos criterios de convergencia"<<endl;
    }
}

void gausSeidel(double *inicial,double mat[maxN][maxN+1], int tam, int iteracoes, double epsilon)
{
    int cont=1,k=0;
    double soma=0;
    double solucao_1[tam];
    for(int i; i<tam; i++)
        solucao_1[i]=0;
    cout<<"{";
    for(int i; i<tam; i++)
    {
        cout<<""<<inicial[i];
        cout<<",";
    }
    cout<<"}";
    cout<<endl;
    cout<<"_______________________________________________"<<endl;

    if(linhas(mat,tam)||sassenfeld(mat,tam))
    {
        condiciona(mat,tam);
        for(int i=0; i<tam; i++)
        {
            soma = mat[i][tam];
            for(int j=tam-1; j >= 0; j--)
            {
                if(j!=i)
                {
                    soma = soma - (mat[i][j])*inicial[j];
                }
            }
            inicial[i] = soma/(mat[i][i]);
            printf("Solucao(%d): X%d = %.5f   ",k,i+1,inicial[i]);
            if(cont == 3)
            {
                double err = erro(solucao_1, inicial,tam);
                printf("Erro: %.5f  ",err);
                cout<<endl;
                cont=0;
            }
            cont++;
        }
        k++;

        do
        {
            for(int i=0; i<tam; i++)
                solucao_1[i]=inicial[i];

            for(int i=0; i<tam; i++)
            {
                soma = mat[i][tam];
                for(int j=tam-1; j >= 0; j--)
                {
                    if(j!=i)
                    {
                        soma = soma - (mat[i][j])*inicial[j];
                    }
                }

                inicial[i] = soma/(mat[i][i]);
                printf("Solucao(%d): X%d = %.5f   ",k,i+1,inicial[i]);

                if(cont == 3)
                {
                    double err = erro(solucao_1, inicial,tam);
                    printf("Erro: %.5f  ",err);
                    cout<<endl;
                    cont=0;
                }
                cont++;
            }
            k++;

        }
        while( erro(solucao_1, inicial,tam)>epsilon && k<iteracoes);
        cout<<endl;
        printf("\nSolucoes:\n");
        for(int i=0; i<tam; i++)
            printf("Solucao %d = %.5f ",i+1,inicial[i]);
        cout<<endl;
    }
    else
    {
        cout<<"Não atende aos criterios de convergencia"<<endl;
    }
}
void eliminacaoGauss(double mat[maxN][maxN+1], int tam)
{
    long double m,pivo,aux;
    //int i,j,k;
    for(int i=0; i<tam-1; i++)
    {
        pivo = mat[i][i];
        for(int k=i; k<tam-1; k++)
        {
            if(fabs(mat[k+1][i]) > pivo)
            {
                printf("troca da linha %d com a linha %d \n",k+1,i);
                for(int j=0; j<=tam; j++)
                {
                    aux = mat[k+1][j];
                    mat[k+1][j] = mat[i][j];
                    mat[i][j] = aux;
                }
                for( int l=0; l<tam; l++)
                {
                    for( int m=0; m<tam+1; m++)
                    {
                        printf("%7g| ", mat[l][m]);
                    }
                    cout<<endl;
                }
                cout<<"______________________________________________"<<endl;
                pivo = mat[i][i];
            }
        }

        for(int k=i+1; k<tam; k++)
        {
            m = (mat[k][i])/pivo;
            for(int j=0; j<=tam; j++)
            {
                mat[k][j] = mat[k][j] - m*(mat[i][j]);
            }
        }
        cout<<"Matriz escalonada:"<<endl;
        for( int i=0; i<tam; i++)
        {
            for( int j=0; j<tam+1; j++)
            {
                printf("%7g| ", mat[i][j]);
            }
            cout<<endl;
        }
        cout<<"______________________________________________"<<endl;
    }
    double solucao[tam],soma;
    solucao[tam-1] = (mat[tam-1][tam])/(mat[tam-1][tam-1]);
    int k = 3;
    for(int i=tam-2; i>=0; i--)
    {
        soma = mat[i][tam];
        for(int j=tam-1; j>i; j--)
            soma = soma - (mat[i][j])*(solucao[j]);
        solucao[(tam+1)-k] = soma/(mat[i][i]);
        k++;
    }
    cout << "\nSolucoes:" << endl;
    for(int i=0; i<tam; i++)
        printf(" X%d = %g ",i+1,solucao[i]);
    cout << "\n_______________________________________________________________" << endl;

}
void menu(int funcao)
{
    double mat[maxN][maxN+1];
    int tam;
    switch(funcao)
    {
    case 0:
    {
        tam = lerEliminacaoGauss(mat);
        eliminacaoGauss(mat,tam);
        break;
    }

    case 1:
    {
        tam = lerJacobi(mat);
        double inicial[tam] = {0,0,0};
        jacobi(inicial,mat,tam,100,pow(10,-5));
        break;
    }

    case 2:
    {
        tam = lerSeidel(mat);
        double inicial2[tam] = {0,0,0};
        gausSeidel(inicial2,mat,tam,100,pow(10,-5));
        break;
    }


    }
}
int main()
{
    int eliGaus = 0, gJacobi = 1, gSeidel = 2;
    menu(eliGaus);
    menu(gJacobi);
    menu(gSeidel);

    return 0;
}
