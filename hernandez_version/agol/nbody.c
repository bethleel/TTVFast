#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define PI 3.14159265358979323846
#define GNEWT 39.4845
#define YEAR 365.242
#define loop(idx,last) for (idx=0; idx<last; idx++)
#define NMAX 15
#define NDIM 3
#include "universal.h"

int main(argc,argv)
int argc;
char **argv;
{
    double m[NMAX],x[NDIM][NMAX],v[NDIM][NMAX];
    double E0,L0[NDIM],p0[NDIM],xcm0[NDIM];
    double E,L[NDIM],p[NDIM],xcm[NDIM];
    double t=0,h,tmax;
    FILE *fcons;
    int i=0,n,pair[NMAX][NMAX];
    void infig1(),phi2(),consq();

    /*Choose an output file*/
    fcons = fopen("data/fcons.txt","w");
    h = atof(argv[1]);
    tmax = atof(argv[2]);
    infig1(m,x,v,&n,pair);
    consq(m,x,v,n,&E0,L0,p0,xcm0);
    t = 0;
    while(t < tmax){
      t += h;
      i++;
      phi2(x,v,h,m,n,pair);
      loop(i,n-1){
        fprintf(fcons," %16.12f %16.12f %16.12f %16.12f",x[0][i],x[1][i],v[0][i],v[1][i]);
      }
      fprintf(fcons," %16.12f %16.12f %16.12f %16.12f\n",x[0][n-1],x[1][n-1],v[0][n-1],v[1][n-1]);
    }
    consq(m,x,v,n,&E,L,p,xcm);
    printf("dE/E=%g\n",(E0-E)/E0);
    fclose(fcons);
    exit(0);
}

double sgn(x)
     double x;
{
  if(x>0.0)
    return(1.0);
  else
    return(-1.0);
}

void centerm(m,x,v,delx,delv,i,j,h)
    double m[NMAX],x[NDIM][NMAX],v[NDIM][NMAX],delx[NDIM],delv[NDIM], h;
    int i,j;
{
    double mij,vcm[NDIM];
    int k;

    mij =m[i] + m[j];
    if(mij == 0){
        loop(k,NDIM){
          vcm[k] = (v[k][i]+v[k][j])/h;
          x[k][i] += h*vcm[k];
          v[k][j] += h*vcm[k];
        }
    }else{
      loop(k,NDIM){
        vcm[k] = (m[i]*v[k][i] + m[j]*v[k][j])/mij;
        x[k][i] += m[j]/mij*delx[k] + h*vcm[k];
        x[k][j] += -m[i]/mij*delx[k] + h*vcm[k];
        v[k][i] += m[j]/mij*delv[k];
        v[k][j] += -m[i]/mij*delv[k];
      }
    }
}


void driftij(x,v,i,j,h)
double x[NDIM][NMAX],v[NDIM][NMAX],h;
int i,j;
{
  int k;
  loop(k,NDIM){
    x[k][i] += h*v[k][i];
    x[k][j] += h*v[k][j];
  }
}



void keplerij(m,x,v,i,j,h)
double m[NMAX],x[NDIM][NMAX],v[NDIM][NMAX],h;
int i,j;
{
  double gm,delx[NDIM],delv[NDIM];
  int k,keplerm();
  void kepler_step();
  State s0,s;

  s0.x = x[0][i] - x[0][j];
  s0.y = x[1][i] - x[1][j];
  s0.z = x[2][i] - x[2][j];
  s0.xd = v[0][i] - v[0][j];
  s0.yd = v[1][i] - v[1][j];
  s0.zd = v[2][i] - v[2][j];
  gm = GNEWT*(m[i]+m[j]);
  if(gm == 0){
    loop(k,NDIM){
      x[k][i] += h*v[k][i];
      x[k][j] += h*v[k][j];
    }
  }else{
    /*different motions if at least one particle massive*/
    kepler_step(gm, h, &s0, &s);
    delx[0] = s.x - s0.x;
    delx[1] = s.y - s0.y;
    delx[2] = s.z - s0.z;
    delv[0] = s.xd - s0.xd;
    delv[1] = s.yd - s0.yd;
    delv[2] = s.zd - s0.zd; 
    /*go back absolute coord, add cm motion*/
    centerm(m,x,v,delx,delv,i,j,h);
  }
}


void drift(x,v,h,n)
double x[NDIM][NMAX],v[NDIM][NMAX],h;
int n;
{
  int i,j;

  loop(i,n){
    loop(j,NDIM){
      x[j][i] += h*v[j][i];
    }
  }
}



void kickfast(x,v,h,m,n,pair)
double x[NDIM][NMAX],v[NDIM][NMAX],h,m[NMAX];
int n,pair[NMAX][NMAX];
{
  double rij[NDIM],r2,r3,fac;
  int i,j,k;

  loop(i,n){
    for(j=i+1; j<n; j++){
      if(pair[i][j] == 1){
        r2 = 0;
        loop(k,NDIM){
          rij[k] = x[k][i] - x[k][j];
          r2 += rij[k]*rij[k];
        }
        r3 = r2*sqrt(r2);
        loop(k,NDIM){ 
          fac = h*GNEWT*rij[k]/r3;
          v[k][i] -= m[j]*fac;
          v[k][j] += m[i]*fac;
        } 
      }
    }
  }
}

void phi2(x,v,h,m,n,pair)
double x[NDIM][NMAX],v[NDIM][NMAX],h,m[NMAX];
int n,pair[NMAX][NMAX];
{
  int i,j;

  drift(x,v,h/2,n);
  kickfast(x,v,h/2,m,n,pair);
  for(i=0; i<n-1; i++){
    for(j=i+1; j<n; j++){
      if(pair[i][j] != 1){
        driftij(x,v,i,j,-h/2);
        keplerij(m,x,v,i,j,h/2);
      }
    }
  }
  for(i=n-2; i>=0;i--){
    for(j=n-1; j>i;j--){
      if(pair[i][j] != 1){
        keplerij(m,x,v,i,j,h/2);
        driftij(x,v,i,j,-h/2);
      }
    }
  }
  kickfast(x,v,h/2,m,n,pair);
  drift(x,v,h/2,n);
}




void infig1(m,x,v,n,pair)
double m[NMAX],x[NDIM][NMAX],v[NDIM][NMAX];
int *n,pair[NMAX][NMAX];
{
  /*A good time step for this problem is 0.001 */
  double a0,a1,e0,e1,mu0,mt0,mu1,mt1,x0,x1,v0,v1,fa0,fb0,fa1,fb1;
  int i,j;


  *n = 3;
  loop(i,*n){
    for(j = i+1;j<*n; j++){
      pair[i][j] = 0; /*all Kepler */
    }
  }
  a0 = 0.0125;
  a1 = 1.0;
  e0 = 0.6;
  e1 = 0;
  m[0] = 1e-3;
  m[1] = 1e-3;
  m[2] = 1.0;
  mu0 = (m[0]*m[1])/(m[0]+m[1]);
  mt0 = m[0]+m[1];
  mu1 = (m[2]*mt0)/(m[2]+mt0);
  mt1 = m[2]+mt0;
  x0 = a0*(1+e0);
  x1 = a1*(1+e1);
  v0 = sqrt(GNEWT*mt0/a0*(1-e0)/(1+e0));
  v1 = sqrt(GNEWT*mt1/a1*(1-e1)/(1+e1));
  fa0 = m[0]/mt0;
  fb0 = m[1]/mt0;
  fa1 = mt0/mt1;
  fb1 = m[2]/mt1;
  loop(i,NMAX){
    loop(j,NDIM){
      x[j][i] = 0;
      v[j][i] = 0;
    }
  }
  x[0][0] = -fb1*x1-fb0*x0;
  x[0][1] = -fb1*x1+fa0*x0;
  x[0][2] = fa1*x1;
  v[1][0] = -fb1*v1-fb0*v0;
  v[1][1] = -fb1*v1+fa0*v0;
  v[1][2] = fa1*v1;
}





void consq(m,x,v,n,Ep,L,p,xcm)
double m[NMAX],x[NDIM][NMAX],v[NDIM][NMAX],p[NDIM],xcm[NDIM],L[NDIM];
double *Ep;
int n;
{
  int i,j,k;
  double E=0,rij,mt=0;

  loop(i,NDIM){
    p[i] = 0;
    xcm[i] = 0;
    L[i] = 0;
  }

  loop(i,n){
    if(m[i] != 0){
      loop(j,NDIM){
        p[j] += m[i]*v[j][i];
        xcm[j] += m[i]*x[j][i];
      }
      L[0] += m[i]*(x[1][i]*v[2][i]-x[2][i]*v[1][i]);
      L[1] += m[i]*(x[2][i]*v[0][i]-x[0][i]*v[2][i]);
      L[2] += m[i]*(x[0][i]*v[1][i]-x[1][i]*v[0][i]);
    }
    loop(j,NDIM){
      E += 1.0/2.0*m[i]*v[j][i]*v[j][i];
    }
    for(j=i+1; j<n; j++){
      if(m[j] != 0){ /*to save time, not necessary*/
        rij = 0;
        loop(k,NDIM){
          rij += (x[k][i]-x[k][j])*(x[k][i]-x[k][j]);
        }
        rij = sqrt(rij);
        E -= GNEWT*m[i]*m[j]/rij;
      }
    }
  }
  loop(i,n){ mt += m[i];}
  if(mt != 0){
    loop(i,NDIM){
      xcm[i] /= mt;
      p[i] /= mt;
    }
  }
  *Ep = E;
}



