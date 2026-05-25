/* homoclinic.c — main solver for the fractional ground state equation
 *
 *   D^gamma_|y| u  -  u  +  u^p  =  0,   1 <= gamma <= 2,  p > 1,
 *
 * with u(+/- infinity) = 0, u > 0, u(y) = u(-y). At gamma = 2 the
 * solution is the classical homoclinic orbit (3a); for 1 <= gamma < 2 it
 * is a fractional ground state with an algebraic tail ~ y^-(gamma+1).
 *
 * Method: pseudospectral FFT on the auxiliary Gierer-Meinhardt system
 * (5a,5b), Crank-Nicolson + Adams-Bashforth time stepping. Boundary
 * truncation error is controlled by a tail condition imposed over a
 * sizeable portion of the domain (see §1 of the paper).
 *
 * This file implements the second-order, second-derivative tail scheme
 * (eq. 9c) and is the variant Kyle's experiments settled on. The four
 * sibling files first_order_*.c / second_order_*.c implement schemes
 * (7b), (7c), (9b), (9c) in isolation for comparison.
 *
 * Input:  homoclinic.in        (12 whitespace-delimited values; see
 *                               homoclinic.in.annotated for descriptions)
 * Output: homoclinic.out       (normalised u on the mesh)
 *         ic.in                (final u and v for warm restarts)
 *
 * Reference: Booker & Nec, "On Accuracy of Numerical Solution to
 * Boundary Value Problems on Infinite Domains with Slow Decay",
 * Math. Model. Nat. Phenom. (2019). DOI 10.1051/mmnp/2019057
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define pi acos(-1.0)

double* allocate(double* P,int n)
{
 /* allocate one one-dimensional array */
 if((P=(double*)malloc(n*sizeof(double)))==NULL){
   printf("\n Out of memory");
   exit(1);
 }
 return P;
}

void memory(double* *U,double* *Un,double* *Uf,
 double* *Unf,double* *Unfp,int N)
{
 int N2=2*N;
 *U=allocate(*U,N);
 *Un=allocate(*Un,N);
 *Uf=allocate(*Uf,N2);
 *Unf=allocate(*Unf,N2);
 *Unfp=allocate(*Unfp,N2);
}

void free_memory(double* *U,double* *Un,double* *Uf,
 double* *Unf,double* *Unfp)
{
 free(*U);
 free(*Un);
 free(*Uf);
 free(*Unf);
 free(*Unfp);
}

void fft(double* Y,int nft,int isign)
{
 int i,j,n,istep,m,mmax;
 double wr,wi,wpr,wpi,wtmp,theta,tmpr,tmpi;
 n=2*nft;
 j=1;
 for(i=1;i<=n;i+=2){
   if(j>i){
     tmpr=*(Y+j-1);
     tmpi=*(Y+j);
     *(Y+j-1)=*(Y+i-1);
     *(Y+j)=*(Y+i);
     *(Y+i-1)=tmpr;
     *(Y+i)=tmpi;
   }
   m=n/2;
   while((m>=2)&&(j>m)){
     j-=m;
     m/=2;
   }
   j+=m;
 }
 mmax=2;
 while(n>mmax){
   istep=2*mmax;
   theta=2.0*pi/(isign*mmax);
   wpr=-2.0*pow(sin(0.5*theta),2);
   wpi=sin(theta);
   wr=1.0;
   wi=0.0;
   for(m=1;m<=mmax;m+=2){
     for(i=m;i<=n;i+=istep){
       j=i+mmax;
       tmpr=wr*(*(Y+j-1))-wi*(*(Y+j));
       tmpi=wr*(*(Y+j))+wi*(*(Y+j-1));
       *(Y+j-1)=*(Y+i-1)-tmpr;
       *(Y+j)=*(Y+i)-tmpi;
       *(Y+i-1)+=tmpr;
       *(Y+i)+=tmpi;
     }
     wtmp=wr;
     wr+=wr*wpr-wi*wpi;
     wi+=wi*wpr+wtmp*wpi;
   }
   mmax=istep;
 }
}

int main()
{
 int N,mt,init,n,ndts,ir,i,j,m;
 double g,p,dt,D,eps,tau;
 double alp,alpg,y,wk,wk2,den,denp,area,conv,th,th_lim;
 double *U,*Un,*Uf,*Unf,*Unfp;
 double *V,*Vn,*Vf,*Vnf,*Vnfp;
 FILE *ic,*out,*prm;
 if((prm=fopen("homoclinic.in","r"))==NULL)
   printf("\n Cannot open prm input file.\n");
 fscanf(prm,"%d %d %lf %lf %d %lf %d %d %d %lf %lf %lf",
   &mt,&init,&g,&p,&n,&dt,&ndts,&ir,&N,&D,&eps,&tau);
 fclose(prm);
 if((out=fopen("homoclinic.out","w"))==NULL)
   printf("\n Cannot open output file.\n");
 conv=0.0; 
 th_lim=pow(10.0,-10.0);

 memory(&U,&Un,&Uf,&Unf,&Unfp,N);
 memory(&V,&Vn,&Vf,&Vnf,&Vnfp,N);

// for(g=1.0;g<=2.0;g+=0.01){ // block for a single homoclinic computation
 
 printf("g=%1.2lf, p=%1.1f\n",g,p);
 alpg=pow(1.0/n,g);
 alp=pow(1.0/n,2.0);

 if(init==0){
   for(j=0;j<N;j++){
     y=2.0*pi*n*(j-N/2)/(N*eps);
     *(U+j)=(pow((p+1.0)/2.0,1.0/(p-1.0))*pow(cosh((p-1.0)/2.0*y),
             -2.0/(p-1.0))-2.0/(1.0+pow(y,2.0)))*
            (g-1.0)+2.0/(1.0+pow(y,2.0));
     *(Uf+2*j)=*(U+j);
     *(Uf+2*j+1)=0.0;
   }
   *V=0.5*(pow(*U,p)+pow(*(U+N-1),p))/(N*eps);
   for(j=1;j<N-1;j++)
     *V+=pow(*(U+j),p)/(N*eps);
   for(j=0;j<N;j++){
     *(V+j)=*V;
     *(Vf+2*j)=*V;
     *(Vf+2*j+1)=0.0;
   }
 }
 else{
   if((ic=fopen("ic.in","r"))==NULL)
     printf("\n Cannot open ic input file.\n");
   for(j=0;j<N;j++){
     fscanf(ic,"%lf %lf",U+j,V+j);
     *(Uf+2*j)=*(U+j);
     *(Uf+2*j+1)=0.0;
     *(Vf+2*j)=*(V+j);
     *(Vf+2*j+1)=0.0;
   }
   fclose(ic);
 }
 for(m=1;m<=mt;m++){
   fft(Uf,N,-1);
   fft(Vf,N,-1);
   for(j=0;j<N;j++){
     *(Un+j)=pow(*(U+j),p)/(*(V+j));
     *(Unf+2*j)=*(Un+j);
     *(Unf+2*j+1)=0.0;
     *(Vn+j)=pow(*(U+j),p);
     *(Vnf+2*j)=*(Vn+j);
     *(Vnf+2*j+1)=0.0;
   }
   fft(Unf,N,-1);
   fft(Vnf,N,-1);
   if(m==1)
     for(j=0;j<2*N;j++){
       *(Unfp+j)=*(Unf+j);
       *(Vnfp+j)=*(Vnf+j);
     }
   for(j=0;j<N-1;j+=2){
     wk=j/2.0;
     wk2=alpg*pow(wk*eps,g);
     den=1.0+0.5*dt*(wk2+1.0);
     denp=1.0-0.5*dt*(wk2+1.0);
     *(Uf+j)=((*(Uf+j))*denp+(1.5*(*(Unf+j))-0.5*(*(Unfp+j)))*dt)/den;
     *(Uf+j+1)=0.0; // Fourier transform of even function is real
     *(Unfp+j)=*(Unf+j);
     *(Unfp+j+1)=*(Unf+j+1);
     wk2=alp*pow(wk,2.0);
     den=1.0+0.5*dt*(D*wk2+1.0)/tau;
     denp=1.0-0.5*dt*(D*wk2+1.0)/tau;
     *(Vf+j)=((*(Vf+j))*denp+(1.5*(*(Vnf+j))-0.5*(*(Vnfp+j)))*dt/
              (eps*tau))/den;
     *(Vf+j+1)=0.0;
     *(Vnfp+j)=*(Vnf+j);
     *(Vnfp+j+1)=*(Vnf+j+1);

     i=j+N;
     wk=(j-N)/2.0;
     wk2=alpg*pow(-wk*eps,g);
     den=1.0+0.5*dt*(wk2+1.0);
     denp=1.0-0.5*dt*(wk2+1.0);
     *(Uf+i)=((*(Uf+i))*denp+(1.5*(*(Unf+i))-0.5*(*(Unfp+i)))*dt)/den;
     *(Uf+i+1)=0.0;
     *(Unfp+i)=*(Unf+i);
     *(Unfp+i+1)=*(Unf+i+1);
     wk2=alp*pow(-wk,2.0);
     den=1.0+0.5*dt*(D*wk2+1.0)/tau;
     denp=1.0-0.5*dt*(D*wk2+1.0)/tau;
     *(Vf+i)=((*(Vf+i))*denp+(1.5*(*(Vnf+i))-0.5*(*(Vnfp+i)))*dt/
              (eps*tau))/den;
     *(Vf+i+1)=0.0;
     *(Vnfp+i)=*(Vnf+i);
     *(Vnfp+i+1)=*(Vnf+i+1);
   }
   fft(Uf,N,1);
   fft(Vf,N,1);
   for(j=0;j<N;j++){
     (*(Uf+2*j))/=N;
     *(Uf+2*j+1)=0.0;
     *(U+j)=*(Uf+2*j);
     (*(Vf+2*j))/=N;
     *(Vf+2*j+1)=0.0;
     *(V+j)=*(Vf+2*j);
   }
   // mend non-monotonicity and negativity at the tails
   i=N/2;

   while(i<N-1 && *(U+i)>0 && *(U+i)>*(U+i+1))
     i++;
   if(i>0.75*N)
     i=0.75*N;

   // decay based on asymptotics
   double A = 10.0;
   double B = (g+1.0)*(A-(g+2.0));
   
   
   for(j=i-4;j<N-2;j++){
     if(g==2.0)
       *(U+j+1)=(*(U+j))*(1.0-2*pi*n/N*tanh((p-1)*pi*n*(j-N/2)/N));
     else
       *(U+j+2)=(*(U+j+1))*(2.0 - A/(j-N/2))-(*(U+j))*(1.0 - A/(j-N/2) + B/((j-N/2)*(j-N/2))); // KYLE'S CHANGE first-order second derivative
   }
   for(j=1;j<N/2;j++) // reflexion
     *(U+j)=*(U+N-j);
   if(g==2.0) // first point of left tail
     *U=(*(U+1))*(1.0-2*pi*n/N*tanh((p-1)*pi*n*(N/2-1)/N));
   else
     *U=(*(U+1))*(1.0-(g+1.0)/(N/2-1)); //SHOULD I CHANGE THIS?
   th=fabs(*(U+N/2)-conv)/(*(U+N/2));
   conv=*(U+N/2);
   if((m%ndts==0 || th<th_lim) && ir==1)
     printf("step=%d, U(0)=%5.25lf\n",m,*(U+N/2));
   if (m==mt || th<th_lim){
     if((ic=fopen("ic.in","w"))==NULL)
       printf("\n Cannot open ic file for writing.\n");
     for(j=0;j<N;j++)
       fprintf(ic,"%5.25lf %5.25lf\n",*(U+j),*(V+j));
     fclose(ic);
     area=0.5*(pow(*U,p)+pow(*(U+N-1),p))/(N*eps);
     for(j=1;j<N-1;j++)
       area+=pow(*(U+j),p)/(N*eps);
     area=pow(area,1.0/(p-1.0));
     for(j=0;j<N;j++){
       *(U+j)/=area;
       fprintf(out,"%5.25lf\n",*(U+j));
     }
     break;
   }
 }

// } // end of gamma loop, block for a single homoclinic computation
 
 fclose(out);
 free_memory(&U,&Un,&Uf,&Unf,&Unfp);
 free_memory(&V,&Vn,&Vf,&Vnf,&Vnfp);
 return 0;
}
