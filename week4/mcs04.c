//
// ASTR698:  Monte Carlo Simulations Week 4 (mcs04.c)
//
// Maria Patterson      18Apr08
// last updated: 20Apr08
//
// fitting galaxy spectra with combinations of stellar spectra
// lots of linear algebra- transpose, inverses
// matrix format T[column][row]
// add noise to galaxy spectra again
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
 
#define ARRAYLIM 800
#define FILENMLIM 100
#define LINELIM 100

float T[ARRAYLIM][ARRAYLIM];

int main(){
  FILE *fp;
  int i;
  char MCS1IN[FILENMLIM] = "stellar_spectra.dat";
  char MCS2IN[FILENMLIM] = "ngc3512.dat";
  char MCS3OUT[FILENMLIM] = "ngc3512fit1.dat";
  char newline[LINELIM];
  float wstar[ARRAYLIM],Odwarf[ARRAYLIM],Bdwarf[ARRAYLIM],Adwarf[ARRAYLIM];
  float Fdwarf[ARRAYLIM],Gdwarf[ARRAYLIM],Kdwarf[ARRAYLIM],Mdwarf[ARRAYLIM];
  float Ggiant[ARRAYLIM],Kgiant[ARRAYLIM],Mgiant[ARRAYLIM];
  float wgal[ARRAYLIM],galflux[ARRAYLIM];
  double n;
  int j,k,number;  
  float Tt[ARRAYLIM][ARRAYLIM];
  float TtT[ARRAYLIM][ARRAYLIM];
  float TtF[ARRAYLIM][ARRAYLIM];
  float s[ARRAYLIM];
  float sum, r;
  float G[ARRAYLIM],G2[ARRAYLIM],G4[ARRAYLIM],G6[ARRAYLIM];
  float chisquared;
  float dummy[ARRAYLIM];
  float noisegal[ARRAYLIM];

  //How many templates do you want to use?
  number=3; 


  //READ IN FILES.  GALAXY THEN STELLAR.  
  i=-1;
  if      ((fp = fopen(MCS2IN,"r")) == NULL)      {
    printf("Could not open input file: %s -", MCS2IN);
    abort();
  }
  while (fgets(newline, sizeof(newline), fp))     {
    if(newline[0] != '#'){
      if ((sscanf(newline,"%f %f", &wgal[i],&galflux[i])) == 2){
	i++;
      }
    }
  }
  n=i;
  fclose(fp);
  i=0;
  if      ((fp = fopen(MCS1IN,"r")) == NULL)      {
    printf("Could not open input file: %s -", MCS1IN);
    abort();
  }  
  while (fgets(newline, sizeof(newline), fp))     {
    if(newline[0] != '#'){
      if ((i<n)&&(sscanf(newline,"%f %f %f %f %f %f %f %f %f %f %f",
			 &wstar[i],&Odwarf[i],&Bdwarf[i],&Adwarf[i],&Fdwarf[i],&Gdwarf[i],&Kdwarf[i],&Mdwarf[i],
			 &Ggiant[i],&Kgiant[i],&Mgiant[i])) == 11){
	i++;
      }
    }
  }
  fclose(fp);
  //check to see if read in right
  //printf("%f %f \n%f %f \n%lf\n", wgal[0],wstar[0],wgal[638],wstar[638],n);


  //MAKE T MATRIX USING TEMPLATES.
  for (i=0;i<n;i++){
    T[0][i]=Bdwarf[i];
  }
  for (i=0;i<n;i++){
    T[1][i]=Gdwarf[i];
  }
  for (i=0;i<n;i++){
    T[2][i]=Kdwarf[i];
  }
  /* for (i=0;i<n;i++){ */
/*     T[3][i]=Fdwarf[i]; */
/*   } */
/*   for (i=0;i<n;i++){ */
/*     T[4][i]=Gdwarf[i]; */
/*   } */
/*   for (i=0;i<n;i++){ */
/*     T[5][i]=Kdwarf[i]; */
/*   } */
/*   for (i=0;i<n;i++){ */
/*     T[6][i]=Mdwarf[i]; */
/*   } */
/*   for (i=0;i<n;i++){ */
/*     T[7][i]=Ggiant[i]; */
/*   } */
/*   for (i=0;i<n;i++){ */
/*     T[8][i]=Kgiant[i]; */
/*   } */
/*   for (i=0;i<n;i++){ */
/*     T[9][i]=Mgiant[i]; */
/*   } */
  //check to see if matrix is correct form
  // printf("%f %f %f \n%f %f %f \n", T[0][0],T[1][0],T[2][0],T[0][(int)(n-1)],T[1][(int)(n-1)],T[2][(int)(n-1)]);
  
  //TAKE TRANSPOSE OF T MATRIX.
  for (j=0;j<number;j++){
    for (i=0;i<n;i++){
      Tt[i][j]=T[j][i];
    }
  }
  //check to see if transpose works
  //printf("%f %f \n %f %f \n %f %f \n", Tt[0][0], Tt[(int)(n-1)][0], Tt[0][1], Tt[(int)(n-1)][1],Tt[0][2], Tt[(int)(n-1)][2]);
  
  //MULTIPLY T MATRIX BY ITS TRANSPOSE.
  sum=0;
  for (i=0;i<number;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(Tt[k][j]*T[i][k]);
      }
      TtT[i][j]=sum;
      sum=0;
    }
  }
  //check to see TtT is correct
  //printf("%e %e %e \n%e %e %e \n", TtT[0][0], TtT[1][0],TtT[2][0],TtT[0][2],TtT[1][2], TtT[2][2]); 

  //INVERT TtT MATRIX.
  for (i=0;i<number;i++){
    r=TtT[i][i];
    TtT[i][i]=1;
    for (j=0;j<number;j++){
      TtT[i][j]=(TtT[i][j])/r;
    }
    for (j=0;j<number;j++){
      if (i!=j){
	r=TtT[j][i];
	TtT[j][i]=0;
	for (k=0;k<number;k++){
	  TtT[j][k]=TtT[j][k] - (r*TtT[i][k]);
	}
      }
    }
  }
  //check to see inversion is correct
  //printf("%e %e %e \n%e %e %e \n", TtT[0][0], TtT[1][0],TtT[2][0],TtT[0][2],TtT[1][2], TtT[2][2]); 

  //MULTIPLY INVERSE BY GALAXY SPECTRA
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(Tt[k][j]*galflux[k]);
      }
      TtF[i][j]=sum;
      sum=0;
    }
  }
  //check TtF
  //printf("%e %e %e \n", TtF[0][0], TtF[0][1],TtF[0][2]); 

  //CALCULATE S MATRIX.
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(TtT[k][j]*TtF[i][k]);
      }
      s[j]=sum;
      sum=0;
    }
  }
  //check s matrix
  for(i=0;i<number;i++){
    printf("%f \n", s[i]);
  }


  //OPEN OUTPUT FILE
  if      ((fp = fopen(MCS3OUT,"w")) == NULL)      {
    printf("Could not open output file: %s -", MCS3OUT);
    abort();
  }  
  
  //MAKE FIT TO GALAXY AND PRINT OUT
  sum=0;
 
  for (j=0;j<n;j++){
    for (k=0;k<number;k++){
      sum=sum+(T[k][j]*s[k]);
      // printf("%f \n", sum);
    }
    G[j]=sum;
    sum=0;
    //  fprintf(fp,"%f %f \n",wstar[j], G[j]);
  }

  //GET CHI SQUARED FIT.
  sum=0;
  for (i=0;i<n;i++){
    sum=sum+((G[i]-galflux[i])*(G[i]-galflux[i]));
  }
  chisquared=sum/(639+number);
  printf("no noise chisquared= %f \n", chisquared);
















  //ADD UNIFORMLY DISTRIBUTED NOISE TO GALAXY SPECTRUM
  //-2 to 2 units
  srand(1);
  for (i=0;i<n;i++){
    dummy[i]=((rand()/2147483647.)*4.)-2.;
    noisegal[i]=galflux[i]+dummy[i];
    //printf("%f %f\n",galflux[i],noisegal[i]);
  }
  //MULTIPLY INVERSE BY GALAXY SPECTRA WITH NOISE
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(Tt[k][j]*noisegal[k]);
      }
      TtF[i][j]=sum;
      sum=0;
    }
  }
  //check TtF
  //printf("%e %e %e \n", TtF[0][0], TtF[0][1],TtF[0][2]); 

  //CALCULATE S MATRIX.
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(TtT[k][j]*TtF[i][k]);
      }
      s[j]=sum;
      sum=0;
    }
  }
  //check s matrix
  for(i=0;i<number;i++){
    printf("%f \n", s[i]);
  }
  //MAKE FIT TO GALAXY AND PRINT OUT
  sum=0;
  for (j=0;j<n;j++){
    for (k=0;k<number;k++){
      sum=sum+(T[k][j]*s[k]);
      // printf("%f \n", sum);
    }
    G2[j]=sum;
    sum=0;
  }
  //redo chi squared.
  sum=0;
  chisquared=0;
  for (i=0;i<n;i++){
    sum=sum+((G2[i]-noisegal[i])*(G2[i]-noisegal[i]));
  }
  chisquared=sum/(639+number);
  printf("2 noise chisquared= %f \n", chisquared);
  

  //-4 to 4 units
  srand(1);
  for (i=0;i<n;i++){
    dummy[i]=((rand()/2147483647.)*8.)-4.;
    noisegal[i]=galflux[i]+dummy[i];
    //printf("%f %f\n",galflux[i],noisegal[i]);
  }
  //MULTIPLY INVERSE BY GALAXY SPECTRA WITH NOISE
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(Tt[k][j]*noisegal[k]);
      }
      TtF[i][j]=sum;
      sum=0;
    }
  }
  //check TtF
  //printf("%e %e %e \n", TtF[0][0], TtF[0][1],TtF[0][2]); 

  //CALCULATE S MATRIX.
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(TtT[k][j]*TtF[i][k]);
      }
      s[j]=sum;
      sum=0;
    }
  }
  //check s matrix
  for(i=0;i<number;i++){
    printf("%f \n", s[i]);
  }
  //MAKE FIT TO GALAXY AND PRINT OUT
  sum=0;
  for (j=0;j<n;j++){
    for (k=0;k<number;k++){
      sum=sum+(T[k][j]*s[k]);
      // printf("%f \n", sum);
    }
    G4[j]=sum;
    sum=0;
  }
  //redo chi squared.
  sum=0;
  chisquared=0;
  for (i=0;i<n;i++){
    sum=sum+((G4[i]-noisegal[i])*(G4[i]-noisegal[i]));
  }
  chisquared=sum/(639+number);
  printf("4 noise chisquared= %f \n", chisquared);
  


  //-6 to 6 units
  srand(1);
  for (i=0;i<n;i++){
    dummy[i]=((rand()/2147483647.)*12.)-6.;
    noisegal[i]=galflux[i]+dummy[i];
    //printf("%f %f\n",galflux[i],noisegal[i]);
  }
  //MULTIPLY INVERSE BY GALAXY SPECTRA WITH NOISE
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(Tt[k][j]*noisegal[k]);
      }
      TtF[i][j]=sum;
      sum=0;
    }
  }
  //check TtF
  //printf("%e %e %e \n", TtF[0][0], TtF[0][1],TtF[0][2]); 

  //CALCULATE S MATRIX.
  sum=0;
  for (i=0;i<1;i++){
    for (j=0;j<number;j++){
      for (k=0;k<n;k++){
	sum=sum+(TtT[k][j]*TtF[i][k]);
      }
      s[j]=sum;
      sum=0;
    }
  }
  //check s matrix
  for(i=0;i<number;i++){
    printf("%f \n", s[i]);
  }
  //MAKE FIT TO GALAXY AND PRINT OUT
  sum=0;
  for (j=0;j<n;j++){
    for (k=0;k<number;k++){
      sum=sum+(T[k][j]*s[k]);
      // printf("%f \n", sum);
    }
    G6[j]=sum;
    sum=0;
  }

  //redo chi squared.
  sum=0;
  chisquared=0;
  for (i=0;i<n;i++){
    sum=sum+((G6[i]-noisegal[i])*(G6[i]-noisegal[i]));
  }
  chisquared=sum/(639+number);
  printf("6 noise chisquared= %f \n", chisquared);

  for (i=0;i<n;i++){
    fprintf(fp,"%f %f %f %f %f \n",wstar[i], G[i], G2[i], G4[i], G6[i]);
  }

  return 0;
}
