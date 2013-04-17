/*
Fisica Computacional EDO Runge Kutta cuarto 0rden
Nombre:Jorge Luis Mayorga
Codigo:20111082
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//-----------------------------------------------------//
// Declaracion de funciones auxiliares                 //
//-----------------------------------------------------//
float dy1dx(float x,float y2);
float dy2dx(float x,float y1,float y2);
float *RungeKuttaFourthOrderStep(float x_old,float y1_old,float y2_old,float h);
float *generatevector_empty(float *vector,int n);
void copyvector(float *vectorfrom,float *vectorto,int n);
void generateP(float *P,float *R,int n_points);
void imprimir_en_bash(float *x,float *y1,float *y2,float *y3,int n_points);
void  escribirxy_txt(float *x,float *y1,float *y2,float *y3,int n_points);
void  graficarGNUPLOT();
//-----------------------------------------------------//
// Declaracion de la funcion Main                      //
//-----------------------------------------------------//
int main(int argc, char **argv)
{
 int n_points=0,i=0;
 float *x=NULL,*y1=NULL,*y2=NULL,h=0.0,a=0.0,b=0.0,*R=NULL,*V=NULL,*P=NULL;
 float y1_0=0.0,y2_0=0.0,y1_prima_temp=0.0,y2_prima_temp=0.0;
 float *temp=NULL,n_points_float=0.0;
 float x_old=0.0,y1_old=0.0,y2_old=0.0; 
 //Fijar condiciones inciales donde y1=y(x) y y2=y'(x)//
 h=0.001;
 y1_0=1.0;
 y2_0=0.0;
 a=0.0;
 b=1.0;

 n_points_float=(b-a)/h;
 n_points=(int)n_points_float;

 //Generar vectores vacios y mallocs
 x=generatevector_empty(x,n_points);
 temp=generatevector_empty(temp,3);
 y1=generatevector_empty(y1,n_points);
 y2=generatevector_empty(y2,n_points);
 R=generatevector_empty(R,n_points);
 V=generatevector_empty(V,n_points);
 P=generatevector_empty(P,n_points);
 
 //Inicio condiciones Iniciales a las variables de fase
 y1[0]=y1_0;
 y2[0]=y2_0;
 x[0]=a;
 for(i=1;i<n_points;i++){

   y1_prima_temp=dy1dx(x[i-1],y1[i-1]);
   y2_prima_temp=dy2dx(x[i-1],y1[i-1],y2[i-1]);

   x_old=x[i-1];
   y1_old=y1[i-1];
   y2_old=y2[i-1];

   temp=RungeKuttaFourthOrderStep(x_old,y1_old,y2_old,h);

   x[i]=temp[0];
   y1[i]=temp[1];
   y2[i]=temp[2];
  }
  
  copyvector(y1,R,n_points);
  copyvector(y2,V,n_points);
  generateP(P,R,n_points);
  
  imprimir_en_bash(x,R,V,P,n_points);
  escribirxy_txt(x,R,V,P,n_points);
  graficarGNUPLOT(); 
return 0;
}

//-----------------------------------------------------//
//            FUNCIONES AUXILIARES                     //
//-----------------------------------------------------//
float dy1dx(float x,float y2){
 float retorno=y2;
 return retorno;
}

float dy2dx(float x,float y1,float y2){
 float M=0,m=0,cte=0,gamma=0,pi=0,G=0;
 M=1000;
 cte=0.2142;
 m=0.01;
 gamma=5/3;
 pi=3.1415;
 G=0.305;
 float retorno=(1/m)*(((-G*M*m)/(pow(y1,2)))+(4*pi*cte/(pow(y1,(3*gamma)))));
 return retorno;
}

float *RungeKuttaFourthOrderStep(float x_old,float y1_old,float y2_old,float h){

 float *retorno=NULL;
 float y_prime_1=0.0;
 float y_prime_2=0.0;
 float x_middle=0.0;
 float k1_y1=0.0,k2_y1=0.0,k3_y1=0.0,k4_y1=0.0;
 float k1_y2=0.0,k2_y2=0.0,k3_y2=0.0,k4_y2=0.0;
 float slope=0.0,x1=0.0,x2=0.0,x3=0.0,x4=0.0;
 float y1_step1=0.0,y1_step2=0.0,y1_step3=0.0,slope_y1=0.0;
 float y2_step1=0.0,y2_step2=0.0,y2_step3=0.0,slope_y2=0.0;
 float x_new=0.0,y1_new=0.0,y2_new=0.0;
 
 retorno=generatevector_empty(retorno,3);
 
 k1_y1=dy1dx(x_old,y2_old);
 k1_y2=dy2dx(x_old,y1_old,y2_old);

 //Paso No.1//
 x1=x_old+(h/2.0);
 y1_step1=y1_old + (h/2.0) * k1_y1;
 y2_step1=y2_old + (h/2.0) * k1_y2;

 k2_y1=dy1dx(x1,y2_step1);
 k2_y2=dy2dx(x1,y1_step1,y2_step1);
 
 //Paso No.2//
 x2=x_old+(h/2.0);
 y1_step2=y1_old + (h/2.0) * k2_y1;
 y2_step2=y2_old + (h/2.0) * k2_y2;
  
 k3_y1=dy1dx(x2,y2_step2);
 k3_y2=dy2dx(x2,y1_step2,y2_step2);

 //Paso No.3//
 x3=x_old+(h/1.0);
 y1_step3=y1_old + (h/1.0) * k3_y1;
 y2_step3=y2_old + (h/1.0) * k3_y2;
  
 k4_y1=dy1dx(x3,y2_step3);
 k4_y2=dy2dx(x3,y1_step3,y2_step3);

 //Paso No.4//
 slope_y1=(1.0/6.0)*(k1_y1+2.0*k2_y1+2.0*k3_y1+k4_y1);
 slope_y2=(1.0/6.0)*(k1_y2+2.0*k2_y2+2.0*k3_y2+k4_y2);
 
 x_new=x_old+h;
 y1_new=y1_old+h*slope_y1;
 y2_new=y2_old+h*slope_y2;

 retorno[0]=x_new;
 retorno[1]=y1_new;
 retorno[2]=y2_new;

 return retorno;
}


float *generatevector_empty(float *vector,int n){
 int i=0;

 if(!(vector = malloc(sizeof(float)*n))){
    fprintf(stderr, "Problem with memory allocation");
  }

 for(i=0;i<n;i++){ 
    vector[i]=0.0;
     }

 return vector;
}

void escribirxy_txt(float *x,float *y1,float *y2,float *y3,int n_points){
 FILE *fileout;
 int i=0;
 fileout=fopen("RungeKuttaEdoDataSolution.txt","w");
 for(i=0;i<(n_points);i++){
 fprintf(fileout,"%f %f %f %f\n",x[i],y1[i],y2[i],y3[i]);
 }
}

void imprimir_en_bash(float *x,float *y1,float *y2,float *y3,int n_points){
 int i=0;
 for(i=0;i<(n_points);i++){
    printf(" x[%d]=%f R[%d]=%f V[%d]=%f P[%d]=%f \n",i,x[i],i,y1[i],i,y2[i],i,y3[i]);
  }
}

void copyvector(float *vectorfrom,float *vectorto,int n){
 int i=0;
 for(i=0;i<n;i++){
   vectorto[i]=vectorfrom[i];
  }
}

void generateP(float *P,float *R,int n_points){
 int i=0;
 float cte=0.2142;
 for(i=0;i<n_points;i++){
   P[i]=cte/pow(R[i],5);
  }
}
void  graficarGNUPLOT(){

  FILE *gplot = popen("gnuplot -persist","w");
  fprintf(gplot, "set term png\n");
  fprintf(gplot, "set output 'Plot_GraficaSolucionTarea7.png'\n");
  fprintf(gplot, "set multiplot layout 2,2 rowsfirst \n");
  
  fprintf(gplot, "set title 'R(t)'\n");
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "plot 'RungeKuttaEdoDataSolution.txt' using 1:2\n");
  
  fprintf(gplot, "set title 'V(t)'\n");
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "plot 'RungeKuttaEdoDataSolution.txt' using 1:3\n");
  
  fprintf(gplot, "set title 'P(t)'\n");
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "plot 'RungeKuttaEdoDataSolution.txt' using 1:4\n");
  
  fprintf(gplot, "set title 'R(t),V(t),P(t)'\n");
  fprintf(gplot, "unset key\n");
  fprintf(gplot, "splot 'RungeKuttaEdoDataSolution.txt' using 2:3:4\n");
  fprintf(gplot, "unset multiplot\n");

  close(gplot);
}

