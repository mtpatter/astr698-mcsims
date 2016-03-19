//
// ASTR698:  Monte Carlo Simulations Week 2 (mcs02.c)
//
// Maria Patterson      21Mar08
// last updated: 02Apr08
//
// magnitude cuts
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define ARRAYLIM 2000000
#define FILENMLIM 100
 

void randomnumber(int num, int seed, int a, int b, int m, int r[] );
void XYSUM(double n, float lsqx[], float lsqy[], float *xysum, float *xxsum, float *yysum, float *xsum, float *ysum, float *slope, float *yinter, float *corrcoef, float *stderrslope, float *stderryinter);



FILE *fp;
int a,b,m,seed,num;
//int r[ARRAYLIM],s;
int dummy[ARRAYLIM];
int i,j;
char MCSOUT[FILENMLIM] = "stuff.out";
char MCS2OUT[FILENMLIM] = "morestuff.out";
char MCS3OUT[FILENMLIM] = "mags12.out";   //CHANGE THIS
char MCS4OUT[FILENMLIM] = "mags212.out";  //AND THIS
char MCS5OUT[FILENMLIM]=  "lsq12.out";    //AND THIS
//int decimalval,strsize;
//char stringval[ARRAYLIM];
int number[ARRAYLIM],rad;
float x[ARRAYLIM];
float y[ARRAYLIM];
float z[ARRAYLIM];
float x2[ARRAYLIM];
float y2[ARRAYLIM];
float z2[ARRAYLIM];
float dist[ARRAYLIM];
float ra[ARRAYLIM],dec[ARRAYLIM],r[ARRAYLIM];
float u1[ARRAYLIM],u2[ARRAYLIM];
float k[ARRAYLIM],v0,v1,v2,v3;
float g1[ARRAYLIM],g2[ARRAYLIM];
float appmag[ARRAYLIM],meanlum,dum1,stdlum,meandist,dum2;
float vcirc[ARRAYLIM],absmag2[ARRAYLIM],appmag2[ARRAYLIM];
float lsqx[ARRAYLIM],lsqy[ARRAYLIM];
float xysum,xxsum,yysum,xsum,ysum,slope,yinter,corrcoef,stderrslope,stderryinter;

int main() {
  int count,visible,visible2;
  int magcut;
  float density[2], pi;
  float m1,b1;
  double n;
  
  magcut=12;   //CHANGE THIS
  printf("magnitude limit %i \n", magcut);
  density[0]=0.25;
  rad=100;
  pi=3.14159265;
  num=density[0]*8*(rad*rad*rad);
  count=0;
  visible=0;
  visible2=0;
  dum1=0;
  dum2=0;
  m1=-7.27;
  b1=-20;


  //INPUT
  //CASE 1
  a=25173;
  b=13849;
  m=pow(2,16);
  seed=1;

  /*      Make sure output file exists.
   *      Write out output values.        */
  if ((fp = fopen(MCSOUT, "w")) == NULL) {
    printf("Could not open output file: %s - ", MCSOUT);
    abort();
  }
  fprintf(fp,"# \n# R RA Dec - %s \n# \n", MCSOUT);


  //THROW STUFF IN A CUBE
  //Throw stuff randomly in x from -100 to 100
  //randomnumber(num, seed, a, b, m, r);
  //srand((unsigned)time(NULL));
  srand(1);  
  for(i=0;i<num;i++){
    //x[i]=(float)(((r[i]/65535.)*2*rad)-rad);
    dummy[i]=rand();
    x[i]=(float)(((dummy[i]/2147483647.)*2*rad)-rad);
    if(i==(num-1)){
      //   seed=r[i];
      srand(dummy[i]);   
    }
  }
  //Throw stuff randomly in y from -100 to 100
  //randomnumber(num, seed,a,b,m,r);
  for(i=0;i<num;i++){
    //y[i]=(float)(((r[i]/65535.)*2*rad)-rad);
    dummy[i]=rand();
    y[i]=(float)(((dummy[i]/2147483647.)*2*rad)-rad);
    if(i==(num-1)){
      // seed=r[i];
      srand(dummy[i]); 
    }
  }
  //Throw stuff randomly in z from -100 to 100
  //randomnumber(num, seed,a,b,m,r);
  for(i=0;i<num;i++){
    //z[i]=(float)(((r[i]/65535.)*2*rad)-rad);
    dummy[i]=rand();
    z[i]=(float)(((dummy[i]/2147483647.)*2*rad)-rad);
    if(i==(num-1)){
      // seed=r[i];
      srand(dummy[i]); 
    }
    //fprintf(fp,"%lf %lf %lf \n", x[i],y[i],z[i]);
    
    //MAKE INTO SPHERE
    dist[i]=sqrt((x[i]*x[i])+(y[i]*y[i])+(z[i]*z[i]));
    if (dist[i]<=rad){
      count=count+1;
      // }
      // }
      // printf("%i \n", count);
      //j=0;
      //for(i=0;i<num;i++){
      // if(dist[i]<=100){
      x2[j]=x[i];
      y2[j]=y[i];
      z2[j]=z[i];
      r[j]=sqrt((x[i]*x[i])+(y[i]*y[i])+(z[i]*z[i]));
      //    fprintf(fp,"%f %f %f %i\n",x2[j],y2[j],z2[j],j);
      j=j+1;
    }
  }
  density[1]=(float)count/((4./3.)*pi*(rad*rad*rad));
  printf("density of cube is %f galaxies per cubic mpc \n", density[0]);
  printf("density of sphere is %f galaxies per cubic mpc\ntotal num= %i \n", density[1], num);
  printf("total in sphere= %i\n", count);


  //CONVERT TO RA,DEC,R.
  for(i=0;i<count;i++){
    //radius-already done above
   
    
    ra[i]=((atan2(y2[i],x2[i]))/(pi))*12.+12.;               //RA
    dec[i]=(((asin(z2[i]/r[i]))/(pi/2.))*90.);               //Dec
    fprintf(fp,"%f %f %f \n", r[i],ra[i],dec[i]);
  }


  /*      Make sure output file exists.
   *      Write out output values.        */
  if ((fp = fopen(MCS2OUT, "w")) == NULL) {
    printf("Could not open output file: %s - ", MCS2OUT);
    abort();
  }
  fprintf(fp,"# \n# Apparent Magnitudes - %s \n# \n", MCS2OUT);

 
  //GENERATE LUMINOSITIES
  //Two sets of uniformly distributed randoms between -1 and 1.
  for(i=0;i<count;i++){
    u1[i]=(float)(((rand()/2147483647.)*2)-1);
    u2[i]=(float)(((rand()/2147483647.)*2)-1);
    k[i]=(u1[i]*u1[i])+(u2[i]*u2[i]);
    while(k[i]<0 || k[i]>1){
      u1[i]=(float)(((rand()/2147483647.)*2)-1);
      u2[i]=(float)(((rand()/2147483647.)*2)-1);
      k[i]=(u1[i]*u1[i])+(u2[i]*u2[i]);
    }
    v0=k[i];
    v1=u1[i];
    v2=u2[i];
    v3=sqrt((-2.*log(v0))/v0);
    g1[i]=((v1*v3)*1.5)-20.;       //use for abs mag
    g2[i]=((v2*v3)*.3);            //use for noise
    //CALCULATE APPARENT MAGNITUDES
    //m-M= 5logd -5
    appmag[i]=5.*log10(r[i]*1000000.)-5.+g1[i];
   
    fprintf(fp,"%f \n", appmag[i]);
  }


  //MAKE MAGNITUDE CUTS AND OUTPUT STUFF
  /*      Make sure output file exists.
   *      Write out output values.        */
  if ((fp = fopen(MCS3OUT, "w")) == NULL) {
    printf("Could not open output file: %s - ", MCS3OUT);
    abort();
  }
  fprintf(fp,"# \n# Mag Cuts - %s \n# \n", MCS3OUT);
  fprintf(fp,"# \n# r ra dec absmag appmag \n# \n");
  //Also calculate meanlum
  for(i=0;i<count;i++){
    if(appmag[i]<=magcut){
      fprintf(fp, "%f %f %f %f %f \n", r[i],ra[i],dec[i],g1[i],appmag[i]);
      visible=visible+1;      
      dum1=dum1+g1[i];    //summing up visible luminosities
      dum2=dum2+r[i];     //summing up distances
    }
  }
  meanlum=dum1/visible;
  dum1=0;
  meandist=dum2/visible;
  dum2=0;


  /*      Make sure output file exists.
   *      Write out output values.        */
  if ((fp = fopen(MCS4OUT, "w")) == NULL) {
    printf("Could not open output file: %s - ", MCS4OUT);
    abort();
  }
  fprintf(fp,"# \n# Mag Cuts again - %s \n# \n", MCS4OUT);
  fprintf(fp,"# \n# r ra dec absmag absmag2 appmag \n# \n");
  //Calculate standard deviation and circular velocities.  Revised absolute magnitude with noise.
  //Calculate new apparent magnitudes.
  for(i=0;i<count;i++){
    vcirc[i]= pow(10,((g1[i]-b1)/m1)+2.5);
    absmag2[i]=g1[i]+g2[i];   //adding noise
    appmag2[i]=5.*log10(r[i]*1000000.)-5.+absmag2[i];
    if(appmag[i]<=magcut){
      dum1=dum1+((g1[i]-meanlum)*(g1[i]-meanlum));
      //fprintf(fp, "%f %f %f %f %f %f\n", r[i],ra[i],dec[i],g1[i],absmag2[i],appmag[i]);
    }
  }
  stdlum=sqrt(dum1/visible);
  printf("%i visible galaxies \n", visible);
  printf("mean luminosity is %f magnitudes\n", meanlum);
  printf("standard deviation is %f \n", stdlum);
  printf("mean distance is %f Mpcs\n", meandist);
  printf("Now add noise\n");
  dum1=0;
  meanlum=0;
  meandist=0;
  stdlum=0;

  //OUTPUT SECOND MAGNITUDE CUTS
  for(i=0;i<count;i++){
    if(appmag2[i]<=magcut){
      visible2=visible2+1;
      dum1=dum1+absmag2[i];    //summing up visible luminosities
      dum2=dum2+r[i];     //summing up distances
      fprintf(fp, "%f %f %f %f %f %f\n", r[i],ra[i],dec[i],g1[i],absmag2[i],appmag[i]);
    }
  }
  meanlum=dum1/visible2;
  dum1=0;
  meandist=dum2/visible2;
  dum2=0;
  for(i=0;i<count;i++){
    if(appmag2[i]<=magcut){
      dum1=dum1+((absmag2[i]-meanlum)*(absmag2[i]-meanlum));
    }
  }
  stdlum=sqrt(dum1/visible2);
  printf("%i visible galaxies \n", visible2);
  printf("mean luminosity is %f magnitudes\n", meanlum);
  printf("standard deviation is %f \n", stdlum);
  printf("mean distance is %f Mpcs\n", meandist);




  /*      Make sure output file exists.
   *      Write out output values.        */
  if ((fp = fopen(MCS5OUT, "w")) == NULL) {
    printf("Could not open output file: %s - ", MCS5OUT);
    abort();
  }
  fprintf(fp,"# \n# LSQ input-%s \n# \n", MCS5OUT);
  fprintf(fp,"# \n# x y  \n# \n");
  //NOW DO LSQ FITTING
  i=0;
  for(j=0;j<count;j++){
    if(appmag2[j]<=magcut){
      lsqx[i]=log10(pow(10,((absmag2[j]-b1)/m1)+2.5)) -2.5;
      lsqy[i]=absmag2[j];
     
      fprintf(fp,"%f %f\n", lsqx[i],lsqy[i]);
      i++;
    }
  }

  n=(double)visible2;
  XYSUM(n,lsqx,lsqy,&xysum,&xxsum,&yysum,&xsum,&ysum,&slope,&yinter,&corrcoef,&stderrslope,&stderryinter);
  printf("slope = %f \nyinter = %f \ncorrcoef = %f \nstderrslope = %f \nstderryinter = %f \n", slope, yinter, corrcoef, stderrslope, stderryinter);




  return 0;
 
}

//FUNCTIONS

// random number generator- mixed
void randomnumber(int num, int seed, int a, int b, int m, int r[]){
  int i;
  r[0]= seed;
  for(i=1; i<num; i++){
    r[i]=(a*r[i-1] + b) % m;
  } 

  return;
}

void XYSUM(double n, float lsqx[], float lsqy[], float *xysum, float *xxsum, float *yysum, float *xsum, float *ysum, float *slope, float *yinter, float *corrcoef, float *stderrslope, float *stderryinter)
{
  int i;
  float ssxy;
  float ssxx;
  float ssyy;
  float xbar;
  float ybar;
  float svar;
  /* Calculate the sum of (xy) */
  *xysum = 0;
  *xxsum = 0;
  *yysum = 0;
  *xsum = 0;
  *ysum = 0;
  for (i=0; i<n; i++) {
    *xysum = *xysum + (lsqx[i]*lsqy[i]);
    *xxsum = *xxsum + (lsqx[i]*lsqx[i]);
    *yysum = *yysum + (lsqy[i]*lsqy[i]);
    *xsum = *xsum + lsqx[i];
    *ysum = *ysum + lsqy[i];
  }
  xbar = *xsum/(double)(n);
  ybar = *ysum/(double)(n);
  ssxy = *xysum - ((double)(n)*xbar*ybar);
  ssxx = *xxsum - ((double)(n)*xbar*xbar);
  ssyy = *yysum - ((double)(n)*ybar*ybar);
  *slope = ssxy/ssxx;
  *yinter= ybar - *slope*xbar;
  *corrcoef= (ssxy*ssxy)/(ssxx*ssyy);
  svar = sqrt(   (ssyy - (*slope*ssxy)))/sqrt(((double)(n)-2)   );
  *stderrslope = svar* (1/(sqrt(ssxx)));
  *stderryinter = svar* (sqrt((1/(double)(n))+((xbar*xbar)/(ssxx)) ));
                                                                                                                                                                             
  return;
 
}
