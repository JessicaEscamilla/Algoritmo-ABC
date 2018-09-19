
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <random>
#include <chrono>
#include <graphics.h>
#include <fstream>

using namespace std;
const double M_PI = 3.1415926535897;

struct parametrosABC{
    int SN;
    int Cmax;
    int Clim;
    int idFC;
};

class costFunction{
private:
    double upperLimit;
    double lowerLimit;
    int idFc;
    double bestPositionComputed[2];
    double bestFitComputed;
    double bestFxComputed;
    double bestPosition[2];
    double bestFit;
    double bestFx;
    string CFname;

public:
    costFunction(parametrosABC misParametros);
    void setupperLimit(double newUpperLimit);
    void setlowerLimit(double newLowerLimit);
    void setidFc(int newidFc);

    void setBestFitComputed(double newFit){bestFitComputed=newFit;};
    void setBestComputedPosition(double x1, double x2){bestPositionComputed[0]=x1;bestPositionComputed[1]=x2;};
    void setBestComputedFx(double newFx){bestFxComputed=newFx;};

    double  getBestPositionComputed1(){return bestPositionComputed[0];};
    double  getBestPositionComputed2(){return bestPositionComputed[1];};
    double  getBestFitComputed(){return bestFitComputed;};
    double  getBestFxComputed(){return bestFxComputed;};

    double  getBestPosition1(){return bestPosition[0];};
    double  getBestPosition2(){return bestPosition[1];};
    double  getBestFit(){return bestFit;};
    double  getBestFx(){return bestFx;};

    double  getupperLimit(){return upperLimit;};
    double  getlowerLimit(){return lowerLimit;};
    int     getidFc(){return idFc;};
    string getName(){return CFname;};
    double  computeFx(double x1, double x2, int idFc);
};

class foodSource{
private:
    int cmej;//contador de cantidad de mejoras
    double coord[2];
    double fx;
    double p;
    double fit;

public:
    foodSource();//este constructor inicializa todo en 0
    void randomInitialization(costFunction miFuncionCosto);
    void computeNeighborhood(double xki, int i, costFunction miFuncionCosto);
    void printResume();//imprime en pantalla el resumen de parametros de la clase
    void setP(double sumFit){p=fit/sumFit;};
    void resetCmej(){cmej=0;};
    double getP(){return p;};
    double getCoord1(){return coord[0];};
    double getCoord2(){return coord[1];};
    double getFit(){return fit;};
    double getfx(){return fx;};
    int getCmej(){return cmej;};
};



//ABC conventional functions
void faseInicializacion(foodSource misFuentes[], costFunction miFuncionCosto, int SN);
void faseEmpleadas(foodSource misFuentes[], costFunction miFuncionCosto, int SN);
void faseExploradoras(foodSource misFuentes[], costFunction miFuncionCosto, int SN);
void faseScouts(foodSource misFuentes[], costFunction miFuncionCosto, parametrosABC misParametros);

//miscelaneous functions
parametrosABC loadParameters();
double calculaFit(double fx);
void updateOutputs(costFunction& miFuncionCosto, foodSource misFuentes[], parametrosABC misParametros);//actualiza salida en pantalla y en archivo texto
double myRand(double downLimit, double upLimit);//genera un numero aleatorio con distibucion uniforme en el rango [uppper, lower]
//se debe garantizar que lowerLimit<upperLimit
void ordenaColumnaProbAsc(double vectorPIdF[][2],int SN);//ordena ascendentmente el vector. Toma como referenci ala columna de probabiliaddes
void calculaProbAcc(double vectorPIdF[][2],int SN);//cacula probabilidades acumuladas. Se asume que el vector de entrada esta ordenado ascendentemente
int myRoulette(double value, double vectorPIdF[][2],int SN);//implementa algoritmo de ruleta
void imprimeResumenBusqueda(costFunction miFuncionCosto, double tiempo);
void actualizaFuentesAntiguas(foodSource misFuentes[], foodSource misFuentesAntiguas[], int SN);
//-----grafic functions
void drawFA(foodSource misFuentes[], foodSource misFuentesAntiguas[], int SN, costFunction miFuncionCosto);//dibuja fuentes de alimento
void draw_framework(int refDistance);//dibuja framework
void draw_circle(int y[]);//dibuja un círculo
void myScale(double x[], double upperL, double lowerL, int y[]);//escala del dominio de busqueda a pixeles
void closeGraphicOutputs();
void initGraphicOutputs(costFunction miFuncionCosto);//inicializa salida en pantalla de graficos

//-------------------------------------------------------------
//-------------------------------------------------------------

//--definitions zone---

foodSource::foodSource(){
    cmej=0;//contador de cantidad de mejoras
    coord[0]=0;
    coord[1]=0;
    fx=0;
    p=0;
    fit=0;
}//este constructor inicializa todo en 0

void foodSource::randomInitialization(costFunction miFuncionCosto){
    coord[0]=miFuncionCosto.getlowerLimit()+(myRand(0.0,1.0)*(miFuncionCosto.getupperLimit()-miFuncionCosto.getlowerLimit()));
    coord[1]=miFuncionCosto.getlowerLimit()+(myRand(0.0,1.0)*(miFuncionCosto.getupperLimit()-miFuncionCosto.getlowerLimit()));
    //actualice propiedades derivadas del objeto: fx, fit
    fx=miFuncionCosto.computeFx(coord[0],coord[1], miFuncionCosto.getidFc());
    fit=calculaFit(fx);


}

void foodSource::computeNeighborhood(double xki, int i, costFunction miFuncionCosto){
    double vmi(0), newFx(0), newFit(0);
    double a=1.0;
    //calcular vm en la dimension i
    vmi=coord[i]+(myRand(-a,a)*(coord[i]-xki));
    if(vmi>miFuncionCosto.getupperLimit())
        vmi=miFuncionCosto.getupperLimit();
    else
    {
            if (vmi<miFuncionCosto.getlowerLimit())
            vmi=miFuncionCosto.getlowerLimit();
    }
    //calcula el fitness del nuevo vm
        //calcula fx
        if(i==0)
            newFx=miFuncionCosto.computeFx(vmi,coord[1], miFuncionCosto.getidFc());
        else
            newFx=miFuncionCosto.computeFx(coord[0],vmi, miFuncionCosto.getidFc());

        //calcula fitness
        newFit=calculaFit(newFx);

    //compara fitness(vm) vs fitness(xm) y almacena las mejores coordenadas
    //actualiza contl

    if(newFit>fit)
    {    //actualiza
      fit=newFit;
      fx=newFx;
      cmej=0;
      if(i==0)
            coord[0]=vmi;
      else
            coord[1]=vmi;


      //evaluar si es necesrio actualizar pi en esta linea
    }
    else
        cmej++;
}//implementa ecuacion (6)


 void foodSource::printResume(){


    }//imprime en pantalla el resumen de parametros de la clase


costFunction::costFunction(parametrosABC misParametros){

    idFc=misParametros.idFC;
    bestPositionComputed[0]=10;
    bestPositionComputed[1]=10;
    bestFitComputed=-1;
    switch(misParametros.idFC)
    {
    case 1://Goldstein–Price function
                bestPosition[0]=0;//the best position obtained from wikipedia
                bestPosition[1]=-1;
                bestFx=3;
                upperLimit=2.0;
                lowerLimit=-2.0;
                CFname="Goldstein–Price";
            break;
    case 2://Sphere function
                bestPosition[0]=0;//the best position obtained from wikipedia
                bestPosition[1]=0;
                bestFx=0;
                upperLimit=5.0;
                lowerLimit=-5.0;
                CFname="Sphere";
            break;
    case 3://Eggholder function
                bestPosition[0]=512;//the best position obtained from wikipedia
                bestPosition[1]=404.2319;
                bestFx=-959.6407;
                upperLimit=512.0;
                lowerLimit=-512.0;
                CFname="Eggholder";
            break;
    case 4://Damavandi function
                bestPosition[0]=2;//the best position obtained from wikipedia
                bestPosition[1]=2;
                bestFx=0;
                upperLimit=14;
                lowerLimit=-14;
                CFname="Damavandi";
            break;
    case 5://Ackleys function
                bestPosition[0]=0;//the best position obtained from wikipedia
                bestPosition[1]=0;
                bestFx=0;
                upperLimit=5.0;
                lowerLimit=-5.0;
                CFname="Ackleys";
            break;
    case 6://beales function
                bestPosition[0]=3.0;//the best position obtained from wikipedia
                bestPosition[1]=0.5;
                bestFx=0;
                upperLimit=4.5;
                lowerLimit=-4.5;
                CFname="Beale's";
            break;
    case 7://booth's function
                bestPosition[0]=1;//the best position obtained from wikipedia
                bestPosition[1]=3;
                bestFx=0;
                upperLimit=10.0;
                lowerLimit=-10.0;
                CFname="Booth's";
            break;
    case 8://Matyas function
                bestPosition[0]=0;//the best position obtained from wikipedia
                bestPosition[1]=0;
                bestFx=0;
                upperLimit=10.0;
                lowerLimit=-10.0;
                CFname="Matyas";
            break;
    case 9://levi function
                bestPosition[0]=1;//the best position obtained from wikipedia
                bestPosition[1]=1;
                bestFx=0;
                upperLimit=10.0;
                lowerLimit=-10.0;
                CFname="Levi";
            break;
    case 10://three-hump camel function
                bestPosition[0]=0;//the best position obtained from wikipedia
                bestPosition[1]=0;
                bestFx=0;
                upperLimit=5.0;
                lowerLimit=-5.0;
                CFname="Three-hump camel";
            break;
    case 11://easom function
                bestPosition[0]=M_PI;//the best position obtained from wikipedia
                bestPosition[1]=M_PI;
                bestFx=-1;
                upperLimit=100.0;
                lowerLimit=-100.0;
                CFname="Easom";
            break;
    case 12://schaffer 2 function
                bestPosition[0]=0;//the best position obtained from wikipedia
                bestPosition[1]=0;
                bestFx=0;
                upperLimit=100.0;
                lowerLimit=-100.0;
                CFname="schaffer 2";
            break;
    case 13://schaffer 4 function
                bestPosition[0]=0;//the best position obtained from wikipedia
                bestPosition[1]=1.25313;
                bestFx=0.292579;
                upperLimit=100.0;
                lowerLimit=-100.0;
                CFname="schaffer 4";
            break;

    default:
        cout<<"error in cost-function identifier"<<endl;

    }
    bestFit=calculaFit(bestFx);
}


double costFunction::computeFx(double x1, double x2, int idFc){

    double y, fact1a, fact1b, fact1, fact2a, fact2b, fact2;
    double term1, term2;
    double term1a, term1a1, term1a1num, term1a1den;
    double exp1, expr1, exp2, expr2;
    double sum1, sum2, sum3;
    double booth1,booth2;
    double lev1, lev2f1,lev2f2, lev3f1, lev3f2, lev2, lev3;
    double easf1, easf2, easexp;
    double scnum, scden;
    double sc4num1, sc4num2, sc4den;


    switch(idFc)
    {
    case 1://Goldstein–Price function
            fact1a = pow((x1 + x2 + 1),2);
            fact1b = 19 - 14*x1 + 3*pow(x1,2) - 14*x2 + 6*x1*x2 + 3*pow(x2,2);
            fact1 = 1 + fact1a*fact1b;
            fact2a = pow(2*x1 - 3*x2,2);
            fact2b = 18 - 32*x1 + 12*pow(x1,2) + 48*x2 - 36*x1*x2 + 27*pow(x2,2);
            fact2 = 30 + fact2a*fact2b;
            y = fact1*fact2;
            break;
    case 2://Sphere function
            y= (pow(x1,2)+pow(x2,2));//esfera
            break;
    case 3:// function
            term1 = -(x2+47) * sin(sqrt(abs(x2+x1/2+47)));
            term2 = -x1 * sin(sqrt(abs(x1-(x2+47))));
            y = term1 + term2;
            break;
    case 4://Damavandi function
            term1a1num=sin(M_PI*(x1-2))*sin(M_PI*(x2-2));
            term1a1den=pow(M_PI,2)*(x1-2)*(x2-2);
            term1a1=term1a1num/term1a1den;
            term1a=pow(fabs(term1a1),5);
            term1=1-term1a;
            term2=2+pow(x1-7,2)+(2*pow(x2-7,2));
            y=term1*term2;
            break;
    case 5://Ackleys function
            exp1=(-0.2*sqrt(0.5*(pow(x1,2)+pow(x2,2))));
            expr1=(-20*exp(exp1));
            exp2=0.5*(cos(2*M_PI*x1)+cos(2*M_PI*x2));
            expr2=exp(exp2);
            y=expr1-expr2+2.7182818+20;
            break;
    case 6://Beale's function
            sum1=pow(1.5-x1+x1*x2,2);
            sum2=pow(2.25-x1+x1*pow(y,2),2);
            sum3=pow(2.625-x1+x1*pow(y,3),2);
            y=sum1+sum2+sum3;
            break;
    case 7://Booth's function
            booth1=pow((x1+2*x2-7),2);
            booth2=pow((2*x1+x2-5),2);
            y=booth1+booth2;
            break;
    case 8://matyas function
            y=(0.26*(pow(x1,2)+pow(x2,2)))-(0.48*x1*x2);
            break;
    case 9://levi function
            lev1=(1/2-1/2*cos(2*3*M_PI*x1));
            lev2f1=pow(x1-1,2);
            lev2f2=1+(1/2-1/2*cos(2*3*M_PI*x2));
            lev2= lev2f1*lev2f2;
            lev3f1=pow(x2-1,2);
            lev3f2=1+(1/2-1/2*cos(2*2*M_PI*x2));
            lev3=lev3f1*lev3f2;
            y=lev1+lev2+lev3;
            break;
    case 10://three-hump camel function
            y=2*pow(x1,2)-1.05*pow(x1,4)+pow(x1,6)/6+x1*x2+pow(x2,2);
            break;
    case 11://Easom function
            easf1=-cos(x1)*cos(x2);
            easexp=pow(x1-M_PI,2)+pow(x2-M_PI,2);
            easf2=exp(-easexp);
            y=easf1*easf2;
            break;
    case 12://Schaffer 2 function
            scnum=(1/2-1/2*cos(2*(pow(x1,2)-pow(x2,2))))-0.5;
            scden=pow(1+0.001*(pow(x1,2)+pow(x2,2)),2);
            y=0.5+scnum/scden;
            break;
    case 13://Schaffer 4 function
            sc4num1=sin(abs(pow(x1,2)-pow(x2,2)));
            sc4num2=(0.5+0.5*cos(2.0*sc4num1))-0.5;
            sc4den=pow(1+0.001*(pow(x1,2)+pow(x2,2)),2);
            y=0.5+sc4num2/sc4den;
            break;

    default:
        cout<<"error in cost-function identifier"<<endl;

    }


return y;

}

void faseInicializacion(foodSource misFuentes[], costFunction miFuncionCosto, int SN){
double sumFit(0);
int m(0);
for(m=0;m<SN;m++)
{
    misFuentes[m].randomInitialization(miFuncionCosto);
    sumFit=sumFit+misFuentes[m].getFit();
}
//actualizar pi
for(m=0;m<SN;m++)
    misFuentes[m].setP(sumFit);
}

void faseEmpleadas(foodSource misFuentes[], costFunction miFuncionCosto, int SN){
double sumFit(0), xki(0.0);
int m(0), k(0), i(0);
for(m=0;m<SN;m++)
{
    //select k diferent from m
    do{
        k=rand()%(SN+1);
    }while(k==m);
    i=rand()%2;

    if(i==0)
            xki=misFuentes[k].getCoord1();
    else
            xki=misFuentes[k].getCoord2();

    misFuentes[m].computeNeighborhood(xki,i,miFuncionCosto);
    sumFit=sumFit+misFuentes[m].getFit();
}
//actualizar pi
for(m=0;m<SN;m++)
    misFuentes[m].setP(sumFit);

}


void faseExploradoras(foodSource misFuentes[], costFunction miFuncionCosto, int SN){
double vectorP[SN][2], xki(0.0),  value(0.0), sumFit(0.0);
int index(0);
int k(0), i(0);
//carga del vector de probabilidades
for (int m=0;m<SN;m++)
{
    vectorP[m][0]=misFuentes[m].getP();
    vectorP[m][1]=static_cast<double>(m);
}
ordenaColumnaProbAsc(vectorP,SN);
calculaProbAcc(vectorP,SN);

for (int m=0;m<SN;m++)
{
    value=myRand(0.0,1.0);
    index=myRoulette(value,vectorP,SN);
    //select k diferent from index
    do{
        k=rand()%(SN+1);
    }while(k==index);
    i=rand()%2;
    if(i==0)
            xki=misFuentes[k].getCoord1();
    else
            xki=misFuentes[k].getCoord2();

    misFuentes[m].computeNeighborhood(xki,i,miFuncionCosto);
    sumFit=sumFit+misFuentes[m].getFit();
}

//actualizar pi
for(int m=0;m<SN;m++)
    misFuentes[m].setP(sumFit);


}

void faseScouts(foodSource misFuentes[], costFunction miFuncionCosto, parametrosABC misParametros){
double sumFit(0.0);

for(int m=0;m<misParametros.SN;m++)
{
    if(misFuentes[m].getCmej()>misParametros.Clim)
    {
        misFuentes[m].randomInitialization(miFuncionCosto);
        misFuentes[m].resetCmej();
    }

    sumFit=misFuentes[m].getFit()+sumFit;
}

//actualizar pi
for(int m=0;m<misParametros.SN;m++)
    misFuentes[m].setP(sumFit);

}

void updateOutputs(costFunction& miFuncionCosto, foodSource misFuentes[], parametrosABC misParametros){
int m(0);
int x[2];
for (m=0;m<misParametros.SN;m++)
{
    if(misFuentes[m].getFit()>miFuncionCosto.getBestFitComputed())
    {
        miFuncionCosto.setBestFitComputed(misFuentes[m].getFit());
        miFuncionCosto.setBestComputedPosition(misFuentes[m].getCoord1(),misFuentes[m].getCoord2());
        miFuncionCosto.setBestComputedFx(misFuentes[m].getfx());
    }
}
}

void actualizaFuentesAntiguas(foodSource misFuentes[], foodSource misFuentesAntiguas[], int SN){
int referen = 10;
for (int m=0;m<SN;m++)
    misFuentesAntiguas[m]=misFuentes[m];
    draw_framework(referen);
}

void imprimeResumenBusqueda(costFunction miFuncionCosto, double tiempo){
    cout.setf(ios::fixed);
    cout.setf(ios::showpoint);
    cout.precision(3);
    cout<<"Resumen de resultados para la funciOn "<<miFuncionCosto.getName()<<endl;
    cout   <<"\t\t"<<"c[0]"<<"\t"<<"c[1]"<<"\t"<<"fitness"<<"\t"<<"fx"<<"\t"<<"tiempo"<<endl;
    cout    <<"Computed"<<"\t"<<miFuncionCosto.getBestPositionComputed1()<<"\t"
            <<miFuncionCosto.getBestPositionComputed2()<<"\t"
            <<miFuncionCosto.getBestFitComputed()<<"\t"
            <<miFuncionCosto.getBestFxComputed()<<"\t"
            <<tiempo<<endl;
    cout    <<"Reference"<<"\t"<<miFuncionCosto.getBestPosition1()<<"\t"
            <<miFuncionCosto.getBestPosition2()<<"\t"
            <<miFuncionCosto.getBestFit()<<"\t"
            <<miFuncionCosto.getBestFx()<<endl;
}
//miscelaneous functions

double calculaFit(double fx){
    if(fx>=0)
        return 1/(1+fx);
    else
        return 1+fabs(fx);
}

parametrosABC loadParameters(){
parametrosABC misParametros;
    char basura;
    ifstream entradatxt, cintxt;
    entradatxt.open ("ABCparameters.txt");
    if (entradatxt.fail())
        cout<<"Error al abrir el archivo"<<endl;

    entradatxt>>basura>>basura
      >>misParametros.SN;
    entradatxt>> basura>>basura>>basura>>basura>>basura
      >> misParametros.Cmax;
    entradatxt>> basura>>basura>>basura>>basura>>basura
      >> misParametros.Clim;

    entradatxt.close();

    cintxt.open ("ABCidFC.txt");
    if (cintxt.fail())
        cout<<"Error al abrir el archivo"<<endl;

    cintxt>>basura>>basura>>basura>>basura>>basura
      >>misParametros.idFC;
    cintxt.close();

    return misParametros;
}

double myRand(double downLimit, double upLimit){
using namespace std::chrono;

  unsigned seed = chrono::system_clock::now().time_since_epoch().count()*rand();
  default_random_engine generator (seed);
  uniform_real_distribution<double> distribution(downLimit,upLimit);
	return distribution(generator);
}

int myRoulette(double value, double vectorPIdF[][2],int SN){
bool notAccepted=true;
int index(0);
while(notAccepted) {
		index = rand()%SN;
		if(value-vectorPIdF[index][0] <= 0)
            notAccepted=false;
	}
	return vectorPIdF[index][1];//for rounding errors
}


void calculaProbAcc(double vectorPIdF[][2],int SN){
double previousP(0.0);
    for(int i=0;i<SN;i++)
    {
        vectorPIdF[i][0]=previousP+vectorPIdF[i][0];
        previousP=vectorPIdF[i][0];
    }
}
//------
void ordenaColumnaProbAsc(double vectorPIdF[][2],int SN){
      int i(0), j(0), flag = 1;
      double temp0(0.0),temp1(0.0);
      int numLength = SN;
      for(i = 1; (i <= numLength) && flag; i++)
     {
          flag = 0;
          for (j=0; j < (numLength -1); j++)
         {
               if (vectorPIdF[j+1][0] < vectorPIdF[j][0])      // ascending order simply changes to <
              {
                    temp0 = vectorPIdF[j][0];             // swap elements
                    vectorPIdF[j][0] = vectorPIdF[j+1][0];
                    vectorPIdF[j+1][0] = temp0;

                    temp1 = vectorPIdF[j][1];             // swap elements
                    vectorPIdF[j][1] = vectorPIdF[j+1][1];
                    vectorPIdF[j+1][1] = temp1;
                    flag = 1;
               }
          }
     }
        //arrays are passed to functions by address; nothing is returned
}//ordena ascendentmente el vector. Toma como referenci a la columna de probabiliaddes

void closeGraphicOutputs(){
           /* clean up graphics*/
   getch();
   closegraph();
}
void initGraphicOutputs(costFunction miFuncionCosto){
    int referen;
    int y[2];
    double x[2]={0,0};
    //int xlimits[2]={0,0};
    double upperL(0), lowerL(0);

    //inicializar gráficos
        int gdriver = DETECT, gmode, errorcode;
       initgraph(&gdriver, &gmode, "");
       errorcode = graphresult();
       if (errorcode != grOk) {  // an error occurred
          cout<<"Graphics error: %s\n", grapherrormsg(errorcode);
          cout<<"Press any key to halt:";
          getch();
          exit(1);               // terminate with an error code
       }
    referen=10;
    draw_framework(referen);
    //draw best reported solution in wikipedia
    upperL=miFuncionCosto.getupperLimit();
    lowerL=miFuncionCosto.getlowerLimit();
    x[0]=miFuncionCosto.getBestPosition1();
    x[1]=miFuncionCosto.getBestPosition2();
    myScale(x, upperL, lowerL,y);
    draw_circle(y);
   return;
}

void drawFA(foodSource misFuentes[], foodSource misFuentesAntiguas[], int SN, costFunction miFuncionCosto)
{
    double x[2]={0,0};
    int y[2];
    //borrar fuentes antiguas
    setcolor(0);
    for (int m=0;m<SN; m++)
    {
        x[0]=misFuentesAntiguas[m].getCoord1();
        x[1]=misFuentesAntiguas[m].getCoord2();
        myScale(x, miFuncionCosto.getupperLimit(), miFuncionCosto.getlowerLimit(), y);
        draw_circle(y);
    }
    //escribir fuentes nuevas
    setcolor(YELLOW);
        for (int m=0;m<SN; m++)
    {
        x[0]=misFuentes[m].getCoord1();
        x[1]=misFuentes[m].getCoord2();
        myScale(x, miFuncionCosto.getupperLimit(), miFuncionCosto.getlowerLimit(),y);
        draw_circle(y);
    }
    setcolor(WHITE);
    x[0]=miFuncionCosto.getBestPosition1();
    x[1]=miFuncionCosto.getBestPosition2();
    myScale(x, miFuncionCosto.getupperLimit(), miFuncionCosto.getlowerLimit(),y);
    draw_circle(y);
    circle(y[0], y[1], 6.0);

}


void draw_framework(int refDistance){
   int left, top, right, bottom;
   int maxx, maxy;
   maxx=getmaxx();
   maxy=getmaxy();

   left = refDistance;
   top = refDistance;
   right = maxx - refDistance;
   bottom = maxy - refDistance;
   /* draw a rectangle */
   rectangle(left,top,right,bottom);

   setlinestyle(DOTTED_LINE, 1, 1);
   line(0+refDistance, maxy/2, maxx-refDistance, maxy/2);//eje x
   line(maxx/2, 0+refDistance, maxx/2, maxy-refDistance);//eje y
   setlinestyle(SOLID_LINE, 1, 1);

}

void draw_circle(int y[]){
    circle(y[0], y[1], 5.0);
}

void myScale(double x[], double upperL, double lowerL, int y[]){
int maxx, maxy, right, bottom, refDistance;
refDistance=10.0;
maxx=getmaxx();
maxy=getmaxy();
right = maxx - 2.0*refDistance;
bottom = maxy - 2.0*refDistance;
y[0]=(right*(1-x[0]/lowerL)/2)+refDistance;
y[1]=(bottom*(1-x[1]/upperL)/2)+refDistance;
}

