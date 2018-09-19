#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <random>
#include <chrono>
#include <graphics.h>
#include "lib.h"

using namespace std;

int main()
{   //declaracion de variables
    clock_t t_ini, t_fin;
    srand (time(NULL)); /* initialize random seed: */
    foodSource misFuentes[100], misFuentesAntiguas[100];
    parametrosABC misParametros;
    misParametros=loadParameters();//carga parametros
    costFunction miFuncionCosto(misParametros);
    int numeroIteraciones(0);
    double tiempo, secs;

    t_ini = clock();


    faseInicializacion(misFuentes, miFuncionCosto, misParametros.SN);
    initGraphicOutputs(miFuncionCosto);
    do{
        faseEmpleadas(misFuentes, miFuncionCosto,misParametros.SN);
        faseExploradoras(misFuentes, miFuncionCosto,misParametros.SN);
        faseScouts(misFuentes, miFuncionCosto,misParametros);
        updateOutputs(miFuncionCosto, misFuentes, misParametros);
        drawFA( misFuentes,  misFuentesAntiguas, misParametros.SN, miFuncionCosto);
        actualizaFuentesAntiguas( misFuentes, misFuentesAntiguas, misParametros.SN);
        numeroIteraciones++;
        delay(100);
    }
    while(numeroIteraciones<misParametros.Cmax);

   // cout    <<"ultima generacion"<<endl;
   // cout    <<"cmej"<<"\t"<<"c[0]"<<"\t"<<"c[1]"
    //        <<"\t"<<"fx"<<"\t"<<"p"<<"\t"<<"fit"<<endl;
    for(int m=0;m<misParametros.SN;m++)
        misFuentes[m].printResume();

    t_fin = clock();
    secs = (double)(t_fin-t_ini)/ CLOCKS_PER_SEC;
    tiempo=secs ;

    imprimeResumenBusqueda(miFuncionCosto, tiempo);
    closeGraphicOutputs();
    return 0;
}
