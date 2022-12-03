#include "STECKMAP/Pierre/POP/sfit.i"
#include "ATON_constants.i"

/* COMPILATION of reaction rates for photo-chemistry of hydrogen
   taken from Maselli et al. 2003 (CRASH)
*/

func gammaH0(T){
  /* DOCUMENT
     collisional ionisation rate from Maselli et al. 2003
     in cm^3/s
  */

  res=5.85e-11;
  res*=sqrt(T);
  res*=1./(1.+sqrt(T/1.e5));
  res*=exp(-157809.1/T);
  return res;
};

func alphaAH(T){
/* DOCUMENT
   case A  recombination rate from Hui & Gnedin 1997
   in cm^3/s
*/
  
  lambda=2.*157807./T;
  res=1.269e-13;
  res*=pow(lambda,1.503);
  res/=pow(1.+(lambda/0.522)^0.47,1.923);
  return res;
};

func recombination_cooling_rate_A_H(T){
  /* DOCUMENT
     case A HII recombinatino cooling rate from Hui & Gnedin 1997
     in erg cm^3 s^-1
  */

  lambda=2.*157807./T;
  res=1.778e-29*pow(lambda, 1.965);
  res/=pow(1.+pow(lambda/0.541,0.502),2.697);
  return res;
};


func alphaBH(T){
  /* DOCUMENT
     case B recombination rate from Hui & Gnedin 1997
     in cm^3/s
     checked against table 1 case B of Ferland et al. 1992
  */

  lambda=2.*157807./T;
  res=2.753e-14;
  res*=pow(lambda,1.5);
  res/=pow(1.+(lambda/2.74)^0.407,2.242);
  return res;
};


func betaH(T){
  /* DOCUMENT
     HI collisional ionisation coefficient from Hui & Gnedin 1997
     in cm^3 s^-1 K^3/2 
  */

  lambda=2.*157807./T;
  res=21.11*pow(T,-3./2.)*exp(-lambda/2.)*pow(lambda,-1.089);
  res/=pow(1.+pow(lambda/0.354,0.874),1.01);
  return res;
};

func ksiH0(T){
  /* DOCUMENT
     collisional ionisation cooling from Maselli et al. 2003
     in erg cm^3 s^-1
  */
  res=1.27e-21*sqrt(T)/(1.+pow(T/1.e5,0.5));
  res*=exp(-157809.1/T);
  return res;
};


func etaH0(T){
  /* DOCUMENT
     recombination cooling for H0 (Maselli et al. 2003)
     case A or B or total ?
     in erg cm^3 s^-1
  */

  res=8.7e-27*sqrt(T)*pow(T/1.e3,-0.2)/(1.+pow(T/1.e6,0.7));
  return res;
};

func psiH0(T){
  /* DOCUMENT
     collisinal exciation cooling for H0 (Maselli et al. 2003)
     in erg cm^3 s^-1
  */

  res=7.5e-19/(1.+pow(T/1.e5,0.5));
  res*=exp(-118348./T);
  return res;
};

func betabremsstrahlung(T){
  /* DOCUMENT
     Bremsstrahlung cooling from Maselli et al. 2003
     in erg cm^3 s^-1
     CAREFUL: we took the densities out of the formula
     so one needs to multiply the result by rho_electrons^2 (in case of pure Hydrogen chemistry
  */
  res=1.42e-27*sqrt(T);
  return res;
};

func coolingrate(T,x){
  /* DOCUMENT
     just the sum of the cooling terms we have implemented
     just multyply by rho^2
     in erg cm^3 s^-1
  */

  res=betabremsstrahlung(T)*x^2+psiH0(T)*(1.-x)^2+ksiH0(T)*(1.-x)^2+etaH0(T)*x^2;
  return res;
};



    
  



if(0){
  // PLOT SOME OF THE COEFFS
  nT=10000;
  T=spanl(1.e2,1.e8,nT);
  ws,1;
  plh,betaH(T),T,type=2;
  //  plh,gammaH0(T),T;
  //  plh,alphaAH(T),T,color="blue";
  //  plh,alphaBH(T),T,color="red";
  logxy,1,1;
  xyleg,"T(K)","cm^3 / s";
  pltitle,"collisional ionisation rate / recombination A/B rates";

  // is it ok that the alphaB>alphaA above 10^6 ?
  // shouldnt it always be smaller ?

 };
