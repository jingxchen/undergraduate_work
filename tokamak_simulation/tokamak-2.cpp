
static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol : use a random exact solution vector\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";


#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
extern PetscScalar Br[199][110][9],Bt[199][110][9],Bp[199][110][9],vr[199][110][8],vt[199][110][9],vp[199][110][9],P[199][110][9],rho[199][110][9];
extern PetscScalar r[199][110][9],theta[199][110][9];

PetscInt nphi = 8,npsi = 111,nthe = 101,n2th = 2*(nthe-1);
PetscScalar dtheta = 0.1,dphi = 0.1,dr = 0.1,dt = 0.1,a=4,gamm=1;

/*
f=rho[i][j][k],P[i][j][k],v,B;t=theta[i][j][k];P[i][j][k]=phi;
fri0=fri-1jk;fr=frijk;fri1=fri+1jk;
frj0=frij-1k;fr=frijk;frj1=frij+1k;
frk0=frijk-1;fr=frijk;frk1=frijk+1;
*/

PetscScalar dF1drho(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s=1/dt+(a+2*r[i][j][k]*cos(theta[i][j][k]))/(r[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])))*vr[i][j][k]+(vr[i+1][j][k]-vr[i-1][j][k])/(2*dr)-sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))*vt[i][j][k]+1/r[i][j][k]*(vt[i][j+1][k]-vt[i][j-1][k])/(2*dtheta)+(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi*(a+cos(theta[i][j][k])));
  return -s;
}
PetscScalar dF1drhoi0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s=-vr[i][j][k]/(2*dr);
  return -s;
}
PetscScalar dF1drhoi1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s=vr[i][j][k]/(2*dr);
  return -s;
}
PetscScalar dF1drhoj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s=-vt[i][j][k]/(2*r[i][j][k]*dtheta);
  return -s;
}
PetscScalar dF1drhoj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s=vt[i][j][k]/(2*r[i][j][k]*dtheta);
  return -s;
}
PetscScalar dF1drhok0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*dphi*2);
  return -s;
}
PetscScalar dF1drhok1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*dphi*2);
  return -s;
}
PetscScalar dF1dvri0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s= -rho[i][j][k]/(2*dr);
  return -s;
}
PetscScalar dF1dvri1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s= rho[i][j][k]/(2*dr);
  return -s;
}
PetscScalar dF1dvtj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -rho[i][j][k]/(2*r[i][j][k]*dtheta);
  return -s;
}
PetscScalar dF1dvtj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = rho[i][j][k]/(2*r[i][j][k]*dtheta);
  return -s;
}
PetscScalar dF1dvpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s= -rho[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*2*dphi);
  return -s;
}
PetscScalar dF1dvpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s= rho[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*2*dphi);
  return -s;
}
PetscScalar dF1dvr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s= (a+r[i][j][k]*cos(theta[i][j][k])*2)/(r[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])))*rho[i][j][k]+(rho[i+1][j][k]-rho[i-1][j][k])/(2*dr);
  return -s;
}
PetscScalar dF1dvt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))*rho[i][j][k]+1/r[i][j][k]*(rho[i][j+1][k]-rho[i][j-1][k])/(2*dtheta);
  return -s;
}
PetscScalar dF1dvp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (rho[i][j][k+1]-rho[i][j][k-1])/((a+r[i][j][k]*cos(theta[i][j][k]))*2*dphi);
  return -s;
}
//********************
PetscScalar dF2dP(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-1/dt-gamm/(r[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])))*(r[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k]))*(vr[i+1][j][k]-vr[i-1][j][k])/(2*dr)+(a+2*r[i][j][k]*cos(theta[i][j][k]))*vr[i][j][k]+(a+r[i][j][k]*cos(theta[i][j][k]))*(vt[i][j+1][k]-vt[i][j-1][k])/(2*dtheta)-r[i][j][k]*sin(theta[i][j][k])*vt[i][j][k]+r[i][j][k]*(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi));
	return s;
}
PetscScalar dF2dvr(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s = -(P[i+1][j][k]-P[i-1][j][k])/(2*dr)-gamm*P[i][j][k]/(r[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])))*(a+2*r[i][j][k]*cos(theta[i][j][k]));
	return s;
}
PetscScalar dF2dvt(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-1/r[i][j][k]*(P[i][j+1][k]-P[i][j-1][k])/(2*dtheta)+gamm*P[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*sin(theta[i][j][k]);
	return s;
}
PetscScalar dF2dvp(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-1/(a+r[i][j][k]*cos(theta[i][j][k]))*(P[i][j][k+1]-P[i][j][k-1])/(2*dphi);
	return s;
}
PetscScalar dF2dPi1(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-vr[i][j][k]/(2*dr);
	return s;
}
PetscScalar dF2dPi0(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=vr[i][j][k]/(2*dr);
	return s;
}
PetscScalar dF2dPj1(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-vt[i][j][k]/(2*r[i][j][k]*dtheta);
	return s;
}
PetscScalar dF2dPj0(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=vt[i][j][k]/(2*r[i][j][k]*dtheta);
	return s;
}
PetscScalar dF2dPk1(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-1/(a+r[i][j][k]*cos(theta[i][j][k]))*vp[i][j][k]/(2*dphi);
	return s;
}
PetscScalar dF2dPk0(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=1/(a+r[i][j][k]*cos(theta[i][j][k]))*vp[i][j][k]/(2*dphi);
	return s;
}
PetscScalar dF2dvri1(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-gamm*P[i][j][k]/(2*dr);
	return s;
}
PetscScalar dF2dvri0(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=gamm*P[i][j][k]/(2*dr);
	return s;
}
PetscScalar dF2dvtj1(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-gamm*P[i][j][k]/(r[i][j][k]*2*dtheta);
	return s;
}
PetscScalar dF2dvtj0(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=gamm*P[i][j][k]/(r[i][j][k]*2*dtheta);
	return s;
}
PetscScalar dF2dvpk1(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=-gamm*P[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*2*dphi);
	return s;
}
PetscScalar dF2dvpk0(PetscInt i,PetscInt j,PetscInt k){
	PetscScalar s;
	s=gamm*P[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*2*dphi);
	return s;
}
//************
PetscScalar dF3_1dvr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/dt-(vr[i+1][j][k]-vr[i-1][j][k])/(2*dr);
  return s;
}
PetscScalar dF3_1dvt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/r[i][j][k]*(vr[i][j+1][k]-vr[i][j-1][k])/(2*dtheta);
  return s;
}
PetscScalar dF3_1dvp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(a+r[i][j][k]*cos(theta[i][j][k]))*(vr[i][j][k+1]-vr[i][j][k-1])/(2*dphi);
  return s;
}
PetscScalar dF3_1drho(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(rho[i][j][k]*rho[i][j][k])*(-(P[i+1][j][k]-P[i-1][j][k])/(2*dr)+Bp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Br[i][j][k+1]-Br[i][j][k-1])/(2*dphi)-Bp[i][j][k]*Bp[i][j][k]*cos(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))-Bp[i][j][k]*(Bp[i+1][j][k]-Bp[i-1][j][k])/(2*dr)-Bt[i][j][k]*Bt[i][j][k]/r[i][j][k]-Bt[i][j][k]*(Bt[i+1][j][k]-Bt[i-1][j][k])/(2*dr)+(Br[i][j+1][k]-Br[i][j-1][k])*Bt[i][j][k]/(2*r[i][j][k]*dtheta));
  return s;
}
PetscScalar dF3_1dBt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -2*Bt[i][j][k]/(rho[i][j][k]*r[i][j][k])-(Bt[i+1][j][k]-Bt[i-1][j][k])/(rho[i][j][k]*2*dr)+(Br[i][j+1][k]-Br[i][j-1][k])/(2*r[i][j][k]*dtheta*rho[i][j][k]);
  return s;
}
PetscScalar dF3_1dBp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (Br[i][j][k+1]-Br[i][j][k-1])/(rho[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k]))*2*dphi)-2*Bp[i][j][k]*cos(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))*rho[i][j][k]-(Bp[i+1][j][k]-Bp[i-1][j][k])/(rho[i][j][k]*2*dr);
  return s;
}
PetscScalar dF3_1dvri1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF3_1dvri0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF3_1dvrj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_1dvrj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_1dvrk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_1dvrk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_1dPi1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(rho[i][j][k]*2*dr);
  return s;
}
PetscScalar dF3_1dPi0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/(rho[i][j][k]*2*dr);
  return s;
}
PetscScalar dF3_1dBrj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_1dBrj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_1dBrk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(rho[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])*2*dphi));
  return s;
}
PetscScalar dF3_1dBrk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(rho[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])*2*dphi));
  return s;
}
PetscScalar dF3_1dBti1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/(2*dr*rho[i][j][k]);
  return s;
}
PetscScalar dF3_1dBti0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/(2*dr*rho[i][j][k]);
  return s;
}
PetscScalar dF3_1dBpi1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(2*dr*rho[i][j][k]);
  return s;
}
PetscScalar dF3_1dBpi0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(2*dr*rho[i][j][k]);
  return s;
}
//***************
PetscScalar dF3_2dvr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(vt[i+1][j][k]-vt[i-1][j][k])/(2*dr);
  return s;
}
PetscScalar dF3_2dvt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/dt-(vt[i][j+1][k]-vt[i][j-1][k])/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dvp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(a+r[i][j][k]*cos(theta[i][j][k]))*(vt[i][j][k+1]-vt[i][j][k-1])/(2*dphi);
  return s;
}
PetscScalar dF3_2drho(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(rho[i][j][k]*rho[i][j][k])*(-(P[i][j+1][k]-P[i][j-1][k])/(2*r[i][j][k]*dtheta)+Bt[i][j][k]*Br[i][j][k]/r[i][j][k]+\
        Br[i][j][k]*(Bt[i+1][j][k]-Bt[i-1][j][k])/(2*dr)-Br[i][j][k]/r[i][j][k]*(Br[i][j+1][k]-Br[i][j-1][k])\
        /(2*dtheta)+(Bp[i][j][k]*Bp[i][j][k])*sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))-\
        Bp[i][j][k]/r[i][j][k]*(Bp[i][j+1][k]-Bp[i][j-1][k])/(2*dtheta)
        +Bp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi));
  return s;
}
PetscScalar dF3_2dBr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/rho[i][j][k]*(Bt[i][j][k]/r[i][j][k]+(Bt[i+1][j][k]-Bt[i-1][j][k])/(2*dr)-\
  (Br[i][j+1][k]-Br[i][j-1][k])/(r[i][j][k]*2*dtheta));
  return s;
}
PetscScalar dF3_2dBt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(rho[i][j][k]*r[i][j][k]);
  return s;
}
PetscScalar dF3_2dBp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/rho[i][j][k]*(2*sin(theta[i][j][k])*Bp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))-(Bp[i][j+1][k]-\
    Bp[i][j-1][k])/(2*r[i][j][k]*dtheta)+(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k]))));
  return s;
}
PetscScalar dF3_2dvti1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vr[i][j][k]/2*(dr);
  return s;
}
PetscScalar dF3_2dvti0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vr[i][j][k]/2*(dr);
  return s;
}
PetscScalar dF3_2dvtj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dvtj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dvtk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_2dvtk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_2dPj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dPj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dBrj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dBrj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dBti1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(rho[i][j][k]*2*dr);
  return s;
}
PetscScalar dF3_2dBti0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/(rho[i][j][k]*2*dr);
  return s;
}
PetscScalar dF3_2dBtk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(2*dphi*rho[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_2dBtk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(2*dphi*rho[i][j][k]*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_2dBpj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_2dBpj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(rho[i][j][k]*r[i][j][k]*2*dtheta);
  return s;
}
//******************
PetscScalar dF3_3dvr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(vp[i+1][j][k]-vp[i-1][j][k])/(2*dr);
  return s;
}
PetscScalar dF3_3dvt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(vp[i][j+1][k]-vp[i][j-1][k])/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_3dvp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/dt-1/(a+r[i][j][k]*cos(theta[i][j][k]))*(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi);
  return s;
}
PetscScalar dF3_3drho(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(rho[i][j][k]*rho[i][j][k])*(-(P[i][j][k+1]-P[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])))-Bp[i][j][k]*Bt[i][j][k]*sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))+Bt[i][j][k]/r[i][j][k]*(Bp[i][j+1][k]-Bp[i][j-1][k])/(2*dtheta)-Bt[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi)-Br[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Br[i][j][k+1]-Br[i][j][k-1])/(2*dphi)+cos(theta[i][j][k])*Bp[i][j][k]*Br[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))+Br[i][j][k]*(Bp[i+1][j][k]-Bp[i-1][j][k])/(2*dr));
  return s;
}
PetscScalar dF3_3dBr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/rho[i][j][k]*(-(Br[i][j][k+1]-Br[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])))+cos(theta[i][j][k])*Bp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))+(Bp[i+1][j][k]-Bp[i-1][j][k])/(2*dr));
  return s;
}
PetscScalar dF3_3dBt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/rho[i][j][k]*(-Bp[i][j][k]*sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))+(Bp[i][j+1][k]-Bp[i][j-1][k])/(2*dtheta*r[i][j][k])-(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k]))));
  return s;
}
PetscScalar dF3_3dBp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/rho[i][j][k]*(-Bt[i][j][k]*sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))+Br[i][j][k]*cos(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_3dvpi1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF3_3dvpi0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF3_3dvpj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_3dvpj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_3dvpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_3dvpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_3dPk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/(rho[i][j][k]*2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_3dPk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = 1/(rho[i][j][k]*2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF3_3dBrk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*rho[i][j][k]*2*dphi);
  return s;
}
PetscScalar dF3_3dBrk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*rho[i][j][k]*2*dphi);
  return s;
}
PetscScalar dF3_3dBtk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*rho[i][j][k]*2*dphi);
  return s;
}
PetscScalar dF3_3dBtk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/((a+r[i][j][k]*cos(theta[i][j][k]))*rho[i][j][k]*2*dphi);
  return s;
}
PetscScalar dF3_3dBpi1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(2*dr*rho[i][j][k]);
  return s;
}
PetscScalar dF3_3dBpi0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/(2*dr*rho[i][j][k]);
  return s;
}
PetscScalar dF3_3dBpj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/(r[i][j][k]*rho[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF3_3dBpj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/(r[i][j][k]*rho[i][j][k]*2*dtheta);
  return s;
}
//********************
PetscScalar dF4_1dBr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/dt-1/r[i][j][k]*(vt[i][j+1][k]-vt[i][j-1][k])/(2*dtheta)+sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))*vt[i][j][k]-(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dBt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (vr[i][j+1][k]-vr[i][j-1][k])/(2*r[i][j][k]*dtheta)-sin(theta[i][j][k])*vr[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]));
  return s;
}
PetscScalar dF4_1dBp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (vr[i][j][k+1]-vr[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dvr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (Bt[i][j+1][k]-Bt[i][j-1][k])/(2*r[i][j][k]*dtheta)-sin(theta[i][j][k])*Bt[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))+(Bp[i][j][k+1]-Bp[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dvt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(Br[i][j+1][k]-Br[i][j-1][k])/(2*r[i][j][k]*dtheta)+sin(theta[i][j][k])*Br[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]));
  return s;
}
PetscScalar dF4_1dvp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(Br[i][j][k+1]-Br[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dBrj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vt[i][j][k]/(2*r[i][j][k]*dtheta);
  return s;
}
PetscScalar dF4_1dBrj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vt[i][j][k]/(2*r[i][j][k]*dtheta);
  return s;
}
PetscScalar dF4_1dBrk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dBrk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dBtj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vr[i][j][k]/(2*dtheta*r[i][j][k]);
  return s;
}
PetscScalar dF4_1dBtj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vr[i][j][k]/(2*dtheta*r[i][j][k]);
  return s;
}
PetscScalar dF4_1dBpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vr[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dBpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vr[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dvrj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF4_1dvrj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/(r[i][j][k]*2*dtheta);
  return s;
}
PetscScalar dF4_1dvrk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dvrk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dvtj1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/(2*dtheta*r[i][j][k]);
  return s;
}
PetscScalar dF4_1dvtj0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(2*dtheta*r[i][j][k]);
  return s;
}
PetscScalar dF4_1dvpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_1dvpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
//********************
PetscScalar dF4_2dBr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (vt[i+1][j][k]-vt[i-1][j][k])/(2*dr)+cos(theta[i][j][k])*vt[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]));
  return s;
}
PetscScalar dF4_2dBt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/dt-(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])))-(vr[i+1][j][k]-vr[i-1][j][k])/(2*dr)-cos(theta[i][j][k])*dr/(a+r[i][j][k]*cos(theta[i][j][k]));
  return s;
}
PetscScalar dF4_2dBp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (vt[i][j][k+1]-vt[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dvr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(Bt[i+1][j][k]-Bt[i-1][j][k])/(2*dr)-cos(theta[i][j][k])*Bt[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]));
  return s;
}
PetscScalar dF4_2dvt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (Bp[i][j][k+1]-Bp[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])))+(Br[i+1][j][k]-Br[i-1][j][k])/(2*dr)+cos(theta[i][j][k])*Br[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]));
  return s;
}
PetscScalar dF4_2dvp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dBri1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vt[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dBri0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vt[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dBti1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dBti0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dBtk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dBtk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dBpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vt[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dBpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vt[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dvri1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dvri0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dvti1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dvti0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_2dvtk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dvtk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dvpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
PetscScalar dF4_2dvpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/(2*dphi*(a+r[i][j][k]*cos(theta[i][j][k])));
  return s;
}
//********************
PetscScalar dF4_3dBr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (vp[i+1][j][k]-vp[i-1][j][k])/(2*dphi)+vp[i][j][k]/r[i][j][k];
  return s;
}
PetscScalar dF4_3dBt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (vp[i][j][k+1]-vp[i][j][k-1])/(2*r[i][j][k]*dphi);
  return s;
}
PetscScalar dF4_3dBp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -1/dt-(vr[i+1][j][k]-vr[i-1][j][k])/(2*dr)-vr[i][j][k]/r[i][j][k]-(vt[i][j][k+1]-vt[i][j][k-1])/(2*r[i][j][k]*dphi);
  return s;
}
PetscScalar dF4_3dvr(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(Bp[i+1][j][k]-Bp[i-1][j][k])/(2*dr)-Bp[i][j][k]/r[i][j][k];
  return s;
}
PetscScalar dF4_3dvt(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -(Bp[i][j][k+1]-Bp[i][j][k-1])/(2*dphi*r[i][j][k]);
  return s;
}
PetscScalar dF4_3dvp(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = (Br[i+1][j][k]-Br[i-1][j][k])/(2*dr)+Br[i][j][k]/r[i][j][k]+(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi*r[i][j][k]);
  return s;
}
PetscScalar dF4_3dBri1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dBri0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dBtk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vp[i][j][k]/(r[i][j][k]*2*dphi);
  return s;
}
PetscScalar dF4_3dBtk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vp[i][j][k]/(r[i][j][k]*2*dphi);
  return s;
}
PetscScalar dF4_3dBpi1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dBpi0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vr[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dBpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -vt[i][j][k]/(2*r[i][j][k]*dphi);
  return s;
}
PetscScalar dF4_3dBpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = vt[i][j][k]/(2*r[i][j][k]*dphi);
  return s;
}
PetscScalar dF4_3dvri1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dvri0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dvtk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bp[i][j][k]/(2*dphi);
  return s;
}
PetscScalar dF4_3dvtk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bp[i][j][k]/(2*dphi);
  return s;
}
PetscScalar dF4_3dvpi1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Br[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dvpi0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Br[i][j][k]/(2*dr);
  return s;
}
PetscScalar dF4_3dvpk1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = Bt[i][j][k]/(2*r[i][j][k]*dphi);
  return s;
}
PetscScalar dF4_3dvpk0(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s;
  s = -Bt[i][j][k]/(2*r[i][j][k]*dphi);
  return s;
}

PetscScalar F1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = ;
  return s;
}
PetscScalar F2(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = 0;
  return s;
}
PetscScalar F3_1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = -vr[i][j][k]/dt-vr[i][j][k]*(vr[i+1][j][k]-vr[i-1][j][k])/(2*dr)-vt[i][j][k]/r[i][j][k]* \
                  (vr[i][j+1][k]-vr[i][j-1][k])/(2*dtheta)-vp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))* \
                  (vr[i][j][k+1]-vr[i][j][k-1])/(2*dphi)-(P[i+1][j][k]-P[i-1][j][k])/(2*dr*rho[i][j][k]) \
                  +1/rho[i][j][k]*(Bp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Br[i][j][k+1]-Br[i][j][k-1]) \
                  /(2*dphi)-Bp[i][j][k]*Bp[i][j][k]*cos(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k])) \
                  -Bp[i][j][k]*(Bp[i+1][j][k]-Bp[i-1][j][k])/(2*dr)-Bt[i][j][k]*Bt[i][j][k]/r[i][j][k]-Bt[i][j][k] \
                  *(Bt[i+1][j][k]-Bt[i-1][j][k])/(2*dr)+Bt[i][j][k]/r[i][j][k]*(Br[i][j+1][k]-Br[i][j-1][k])/(2*dtheta));
  return s;
}
PetscScalar F3_2(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = -vt[i][j][k]/dt-vr[i][j][k]*(vt[i+1][j][k]-vt[i-1][j][k])/(2*dr)-vt[i][j][k]/r[i][j][k]* \
                  (vt[i][j+1][k]-vt[i][j-1][k])/(2*dtheta)-vp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))* \
                  (vt[i][j][k+1]-vt[i][j][k-1])/(2*dphi)-(P[i][j+1][k]-P[i][j-1][k])/(2*dtheta*rho[i][j][k]*r[i][j][k]) \
                  +1/rho[i][j][k]*(Bt[i][j][k]*Br[i][j][k]/r[i][j][k]+Br[i][j][k]*(Bt[i+1][j][k]-Bt[i-1][j][k])/(2*dr) \
                  -Br[i][j][k]/r[i][j][k]*(Br[i][j+1][k]-Br[i][j-1][k])/(2*dtheta)+Bp[i][j][k]*Bp[i][j][k]* \
                  sin(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k]))-Bp[i][j][k]/r[i][j][k]*(Bp[i][j+1][k]-Bp[i][j-1][k]) \
                  /(2*dtheta)+Bp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi)) \
  return s;
}
PetscScalar F3_3(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = -vp[i][j][k]/dt-vr[i][j][k]*(vp[i+1][j][k]-vp[i-1][j][k])/(2*dr)-vt[i][j][k]/r[i][j][k]* \
                  (vp[i][j+1][k]-vp[i][j-1][k])/(2*dtheta)-vp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))* \
                  (vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi)-(P[i][j][k+1]-P[i][j][k-1])/(rho[i][j][k]*(a+r[i][j][k]* \
                  cos(theta[i][j][k]))*2*dphi)+1/rho[i][j][k]*(-Bp[i][j][k]*Bt[i][j][k]*sin(theta[i][j][k])/(a+r[i][j][k] \
                  *cos(theta[i][j][k]))+Bt[i][j][k])/r[i][j][k]*(Bp[i][j+1][k]-Bp[i][j-1][k])/(2*dtheta)-Bt[i][j][k] \
                  /(a+r[i][j][k]*cos(theta[i][j][k]))*(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi)-Br[i][j][k]/(a+r[i][j][k] \
                  *cos(theta[i][j][k]))*(Br[i][j][k+1]-Br[i][j][k-1])/(2*dphi)+cos(theta[i][j][k])*Bp[i][j][k] \
                  *Br[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))+Br[i][j][k]*(Bp[i+1][j][k]-Bp[i-1][j][k])/(2*dr));
  return s;
}
PetscScalar F4_1(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = -Br[i][j][k]/dt-vt[i][j][k]/r[i][j][k]*(Br[i][j+1][k]-Br[i][j-1][k])/(2*dtheta)-Br[i][j][k]/r[i][j][k]* \
                  (vt[i][j+1][k]-vt[i][j-1][k])/(2*dtheta)+vr[i][j][k]/r[i][j][k]*(Bt[i][j+1][k]-Bt[i][j-1][k])/(2*dtheta) \
                  +Bt[i][j][k]/r[i][j][k]*(vr[i][j+1][k]-vr[i][j-1][k])/(2*dtheta)+sin(theta[i][j][k])/(a+r[i][j][k]* \
                  cos(theta[i][j][k]))*(vt[i][j][k]*Br[i][j][k]-vr[i][j][k]*Bt[i][j][k])+vr[i][j][k]/(a+r[i][j][k] \
                  *cos(theta[i][j][k]))*(Bp[i][j][k+1]-Bp[i][j][k-1])/(2*dphi)+Bp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k])) \
                  *(vr[i][j][k+1]-vr[i]j[k-1])/(2*dphi)-vp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Br[i][j][k+1] \
                  -Bp[i][j][k-1])/(2*dphi)-Br[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi);
  return s;
}
PetscScalar F4_2(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = -Bt[i][j][k]/dt-vp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi)- \
                  Bt[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi)+vt[i][j][k] \
                  /(a+r[i][j][k]*cos(theta[i][j][k]))*(Bp[i][j][k+1]-Bp[i][j][k-1])/(2*dphi)+Bp[i][j][k]/(a+r[i][j][k] \
                  *cos(theta[i][j][k]))*(vt[i][j][k+1]-vt[i][j][k-1])/(2*dphi)+vt[i][j][k]/(2*dr)*(Br[i+1][j][k]-Br[i-1][j][k]) \
                  +Br[i][j][k]/(2*dr)*(vt[i+1][j][k]-vt[i-1][j][k])-vr[i][j][k]/(2*dr)*(Bt[i+1][j][k]-Bt[i-1][j][k]) \
                  -Bt[i][j][k]/(2*dr)*(vr[i+1][j][k]-vr[i-1][j][k])+cos(theta[i][j][k])/(a+r[i][j][k]*cos(theta[i][j][k])) \
                  *(vt[i][j][k]*Br[i][j][k]-vr[i][j][k]*Bt[i][j][k]);
  return s;
}
PetscScalar F4_3(PetscInt i,PetscInt j,PetscInt k){
  PetscScalar s = -Bp[i][j][k]/dt-vr[i][j][k]/(2*dr)*(Bp[i+1][j][k]-Bp[i-1][j][k])-Bp[i][j][k]/(2*dr)*(vr[i+1][j][k] \
                  -vr[i-1][j][k])+vp[i][j][k]/(2*dr)*(Br[i+1][j][k]-Br[i-1][j][k])+Br[i][j][k]/(2*dr)*(vp[i+1][j][k] \
                  -vp[i-1][j][k])-1/r[i][j][k]*(vr[i][j][k]*Bp[i][j][k]-vp[i][j][k]*Br[i][j][k])+vp[i][j][k]/r[i][j][k] \
                  *(Bt[i][j][k+1]-Bt[i][j][k-1])/(2*dphi)+Bt[i][j][k]/r[i][j][k]*(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi) \
                  -vp[i][j][k]/(a+r[i][j][k]*cos(theta[i][j][k]))*(Br[i][j][k+1]-Br[i][j][k-1])/(2*dphi)-Br[i][j][k] \
                  /(a+r[i][j][k]*cos(theta[i][j][k]))*(vp[i][j][k+1]-vp[i][j][k-1])/(2*dphi);
  return s;
}


int main(int argc,char **args)
{
  Vec            x,b,u;  /* approx solution, RHS, exact solution */
  Mat            A;        /* linear system matrix */
  KSP            ksp;     /* linear solver context */
  PetscRandom    rctx;     /* random number generator context */
  PetscReal      norm;     /* norm of solution error */
  PetscInt       i,j,k,Ii,J,Istart0,Iend0,n,its,nstep;;
  PetscInt       Istart1,Iend1,Istart2,Iend2,Istart3,Iend3,Istart4,Iend4,Istart5,Iend5,Istart6,Iend6,Istart7,Iend7;
  PetscErrorCode ierr;
  PetscBool      flg = PETSC_FALSE;
  PetscScalar    v,coef;
  PetscMPIInt    rank;
  FILE           *fp;
#if defined(PETSC_USE_LOG)
  PetscLogStage stage;
#endif

  PetscInitialize(&argc,&args,(char*)0,help);
  //ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(i,j,k), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);CHKERRQ(ierr);
  //ierr = MatSeqAIJSetPreallocation(A,5,NULL);CHKERRQ(ierr);
  //ierr = MatSeqSBAIJSetPreallocation(A,1,5,NULL);CHKERRQ(ierr);

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if(rank == 0){
    ierr = MatGetOwnershipRange(A,&Istart0,&Iend0);CHKERRQ(ierr);
  }
  if(rank == 1){
    ierr = MatGetOwnershipRange(A,&Istart1,&Iend1);CHKERRQ(ierr);
  }
  if(rank == 2){
    ierr = MatGetOwnershipRange(A,&Istart2,&Iend2);CHKERRQ(ierr);
  }
  if(rank == 3){
    ierr = MatGetOwnershipRange(A,&Istart3,&Iend3);CHKERRQ(ierr);
  }
  if(rank == 4){
    ierr = MatGetOwnershipRange(A,&Istart4,&Iend4);CHKERRQ(ierr);
  }
  if(rank == 5){
    ierr = MatGetOwnershipRange(A,&Istart5,&Iend5);CHKERRQ(ierr);
  }
  if(rank == 6){
    ierr = MatGetOwnershipRange(A,&Istart6,&Iend6);CHKERRQ(ierr);
  }
  if(rank == 7){
    ierr = MatGetOwnershipRange(A,&Istart7,&Iend7);CHKERRQ(ierr);
  }


//组装矩阵A
  if(rank == 0){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos1 = Istart0 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos2 = Istart0 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos3 = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;
          PetscInt pos4 = Istart0 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;
          PetscInt pos5 = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;
          PetscInt pos6 = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;

          PetscInt pos7 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos8 = Istart5 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos9 = Istart5 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;

          PetscInt pos10 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos11 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;
          PetscInt pos12 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;

          PetscInt pos13 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos14 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;
          PetscInt pos15 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;
          if(j == 1){
            pos3 = Istart0 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
            pos11 = Istart6 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
          }
          if(j == n2th-1){
            pos4 = Istart0 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
            pos12 = Istart6 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
          }
          if(k == 1){
            pos5 = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos14 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
          }
          if(k == nphi){
            pos6 = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos15 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
          }
          if(i == 1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            MatSetValues(A,1,&pos,1,&pos2,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart0 + (i+1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else if(i == n2th-1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            MatSetValues(A,1,&pos,1,&pos1,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart0 + (i-3)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else{
            coef = dF1drho(i,j,k);
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = dF1drhoi0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos1,&coef,INSERT_VALUES);
            coef = dF1drhoi1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos2,&coef,INSERT_VALUES);
            coef = dF1drhoj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos3,&coef,INSERT_VALUES);
            coef = dF1drhoj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos4,&coef,INSERT_VALUES);
            coef = dF1drhok0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos5,&coef,INSERT_VALUES);
            coef = dF1drhok1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos6,&coef,INSERT_VALUES);
            coef = dF1dvr(i,j,k);
            MatSetValues(A,1,&pos,1,&pos7,&coef,INSERT_VALUES);
            coef = dF1dvri0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos8,&coef,INSERT_VALUES);
            coef = dF1dvri1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos9,&coef,INSERT_VALUES);
            coef = dF1dvt(i,j,k);
            MatSetValues(A,1,&pos,1,&pos10,&coef,INSERT_VALUES);
            coef = dF1dvtj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos11,&coef,INSERT_VALUES);
            coef = dF1dvtj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos12,&coef,INSERT_VALUES);
            coef = dF1dvp(i,j,k);
            MatSetValues(A,1,&pos,1,&pos13,&coef,INSERT_VALUES);
            coef = dF1dvpk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos14,&coef,INSERT_VALUES);
            coef = dF1dvpk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos15,&coef,INSERT_VALUES);
          }
        }
      }
    }
    //i==1
    //i==n2th-1
  }
  if(rank == 1){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos1 = Istart1 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos2 = Istart1 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos3 = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;
          PetscInt pos4 = Istart1 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;
          PetscInt pos5 = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;
          PetscInt pos6 = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;

          PetscInt pos7 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos8 = Istart5 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos9 = Istart5 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;

          PetscInt pos10 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos11 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;
          PetscInt pos12 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;

          PetscInt pos13 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          PetscInt pos14 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;
          PetscInt pos15 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;
          if(j == 1){
            pos3 = Istart1 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
            pos11 = Istart6 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
          }
          if(j == n2th-1){
            pos4 = Istart1 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
            pos12 = Istart6 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
          }
          if(k == 1){
            pos5 = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos14 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
          }
          if(k == nphi){
            pos6 = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos15 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
          }
          if(i == 1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            MatSetValues(A,1,&pos,1,&pos2,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart1 + (i+1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else if(i == n2th-1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            MatSetValues(A,1,&pos,1,&pos1,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart1 + (i-3)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else{
            coef = dF2dP(i,j,k);
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = dF2dPi0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos1,&coef,INSERT_VALUES);
            coef = dF2dPi1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos2,&coef,INSERT_VALUES);
            coef = dF2dPj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos3,&coef,INSERT_VALUES);
            coef = dF2dPj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos4,&coef,INSERT_VALUES);
            coef = dF2dPk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos5,&coef,INSERT_VALUES);
            coef = dF2dPk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos6,&coef,INSERT_VALUES);
            coef = dF2dvr(i,j,k);
            MatSetValues(A,1,&pos,1,&pos7,&coef,INSERT_VALUES);
            coef = dF2dvri0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos8,&coef,INSERT_VALUES);
            coef = dF2dvri1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos9,&coef,INSERT_VALUES);
            coef = dF2dvt(i,j,k);
            MatSetValues(A,1,&pos,1,&pos10,&coef,INSERT_VALUES);
            coef = dF2dvtj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos11,&coef,INSERT_VALUES);
            coef = dF2dvtj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos12,&coef,INSERT_VALUES);
            coef = dF2dvp(i,j,k);
            MatSetValues(A,1,&pos,1,&pos13,&coef,INSERT_VALUES);
            coef = dF2dvpk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos14,&coef,INSERT_VALUES);
            coef = dF2dvpk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos15,&coef,INSERT_VALUES);
          }
        }
      }
    }
  }
  if(rank == 2){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Br(i,j,k)
          PetscInt pos1 = Istart2 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;//Br(i,j+1,k)
          PetscInt pos2 = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;//Br(i,j-1,k)
          PetscInt pos3 = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//Br(i,j,k+1)
          PetscInt pos4 = Istart2 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;//Br(i,j,k-1)

          PetscInt pos5 = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Bt(i,j,k)
          PetscInt pos6 = Istart3 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k;//Bt(i,j+1,k)
          PetscInt pos7 = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;//Bt(i,j-1,k)

          PetscInt pos8 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Bp(i,j,k)
          PetscInt pos9 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//Bp(i,j,k+1)
          PetscInt pos10 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;//Bp(i,j,k-1)

          PetscInt pos11 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vr(i,j,k)
          PetscInt pos12 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;//vr(i,j+1,k)
          PetscInt pos13 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;//vr(i,j-1,k)
          PetscInt pos14 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//vr(i,j,k+1)
          PetscInt pos15 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;//vr(i,j,k-1)

          PetscInt pos16 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vt(i,j,k)
          PetscInt pos17 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;//vt(i,j+1,k)
          PetscInt pos18 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-2)*(nphi) + k-1;//vt(i,j-1,k)

          PetscInt pos19 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vp(i,j,k)
          PetscInt pos20 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//vp(i,j,k+1)
          PetscInt pos21 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;//vp(i,j,k-1)
          if(j == 1){
            pos2 = Istart2 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
            pos7 = Istart3 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
            pos13 = Istart5 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
            pos18 = Istart6 + (i-1)*(n2th-1)*(nphi) + (n2th-2)*(nphi) + k-1;
          }
          if(j == n2th-1){
            pos1 = Istart2 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
            pos6 = Istart3 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
            pos12 = Istart5 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
            pos17 = Istart6 + (i-1)*(n2th-1)*(nphi) + (1-1)*(nphi) + k-1;
          }
          if(k == 1){
            pos4 = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos10 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos15 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos21 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
          }
          if(k == nphi){
            pos3 = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos9 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos14 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos20 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
          }
          if(i == 1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            PetscInt posi1 = Istart2 + i*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi1,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart2 + (i+1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else if(i == n2th-1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            PetscInt posi1 = Istart2 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi1,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart2 + (i-3)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else{
            coef = dF4_1dBr(i,j,k);
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = dF4_1dBrj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos1,&coef,INSERT_VALUES);
            coef = dF4_1dBrj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos2,&coef,INSERT_VALUES);
            coef = dF4_1dBrk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos3,&coef,INSERT_VALUES);
            coef = dF4_1dBrk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos4,&coef,INSERT_VALUES);
            coef = dF4_1dBt(i,j,k);
            MatSetValues(A,1,&pos,1,&pos5,&coef,INSERT_VALUES);
            coef = dF4_1dBtj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos6,&coef,INSERT_VALUES);
            coef = dF4_1dBtj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos7,&coef,INSERT_VALUES);
            coef = dF4_1dBp(i,j,k);
            MatSetValues(A,1,&pos,1,&pos8,&coef,INSERT_VALUES);
            coef = dF4_1dBpk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos9,&coef,INSERT_VALUES);
            coef = dF4_1dBpk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos10,&coef,INSERT_VALUES);
            coef = dF4_1dvr(i,j,k);
            MatSetValues(A,1,&pos,1,&pos11,&coef,INSERT_VALUES);
            coef = dF4_1dvrj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos12,&coef,INSERT_VALUES);
            coef = dF4_1dvrj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos13,&coef,INSERT_VALUES);
            coef = dF4_1dvrk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos14,&coef,INSERT_VALUES);
            coef = dF4_1dvrk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos15,&coef,INSERT_VALUES);
            coef = dF4_1dvt(i,j,k);
            MatSetValues(A,1,&pos,1,&pos16,&coef,INSERT_VALUES);
            coef = dF4_1dvtj1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos17,&coef,INSERT_VALUES);
            coef = dF4_1dvtj0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos18,&coef,INSERT_VALUES);
            coef = dF4_1dvp(i,j,k);
            MatSetValues(A,1,&pos,1,&pos19,&coef,INSERT_VALUES);
            coef = dF4_1dvpk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos20,&coef,INSERT_VALUES);
            coef = dF4_1dvpk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos21,&coef,INSERT_VALUES);
          }
        }
      }
    }
  }
  if(rank == 3){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Bt(i,j,k)
          PetscInt pos1 = Istart3 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Bt(i+1,j,k)
          PetscInt pos2 = Istart3 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Bt(i-1,j,k)
          PetscInt pos3 = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//Bt(i,j,k+1)
          PetscInt pos4 = Istart3 + (i-1)*(n2th-1)*(nphi) + (j)*(nphi) + k-1;//Bt(i,j,k-1)

          PetscInt pos5 = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Br(i,j,k)
          PetscInt pos6 = Istart2 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//Br(i+1,j,k)
          PetscInt pos7 = Istart2 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Br(i-1,j,k)

          PetscInt pos8 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//Bp(i,j,k)
          PetscInt pos9 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//Bp(i,j,k+1)
          PetscInt pos10 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;//Bp(i,j,k-1)

          PetscInt pos11 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vt(i,j,k)
          PetscInt pos12 = Istart6 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vt(i+1,j,k)
          PetscInt pos13 = Istart6 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vt(i-1,j,k)
          PetscInt pos14 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//vt(i,j,k+1)
          PetscInt pos15 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;//vt(i,j,k-1)

          PetscInt pos16 = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vr(i,j,k)
          PetscInt pos17 = Istart5 + (i)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vr(i+1,j,k)
          PetscInt pos18 = Istart5 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vr(i-1,j,k)

          PetscInt pos19 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;//vp(i,j,k)
          PetscInt pos20 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k;//vp(i,j,k+1)
          PetscInt pos21 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-2;//vp(i,j,k-1)
          if(k == 1){
            pos4 = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos10 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos15 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
            pos21 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + nphi-1;
          }
          if(k == nphi){
            pos3 = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos9 = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos14 = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
            pos20 = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + 1-1;
          }
          if(i == 1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            PetscInt posi1 = Istart3 + i*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi1,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart3 + (i+1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else if(i == n2th-1){
            coef = 1;
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = 4./3.;
            PetscInt posi1 = Istart3 + (i-2)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi1,&coef,INSERT_VALUES);
            coef = 1./3.;
            PetscInt posi2 = Istart3 + (i-3)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            MatSetValues(A,1,&pos,1,&posi2,&coef,INSERT_VALUES);
          }
          else{
            coef = dF4_2dBt(i,j,k);
            MatSetValues(A,1,&pos,1,&pos,&coef,INSERT_VALUES);
            coef = dF4_2dBti1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos1,&coef,INSERT_VALUES);
            coef = dF4_2dBti0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos2,&coef,INSERT_VALUES);
            coef = dF4_2dBtk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos3,&coef,INSERT_VALUES);
            coef = dF4_2dBtk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos4,&coef,INSERT_VALUES);
            coef = dF4_2dBr(i,j,k);
            MatSetValues(A,1,&pos,1,&pos5,&coef,INSERT_VALUES);
            coef = dF4_2dBri1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos6,&coef,INSERT_VALUES);
            coef = dF4_2dBri0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos7,&coef,INSERT_VALUES);
            coef = dF4_2dBp(i,j,k);
            MatSetValues(A,1,&pos,1,&pos8,&coef,INSERT_VALUES);
            coef = dF4_2dBpk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos9,&coef,INSERT_VALUES);
            coef = dF4_2dBpk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos10,&coef,INSERT_VALUES);
            coef = dF4_2dvt(i,j,k);
            MatSetValues(A,1,&pos,1,&pos11,&coef,INSERT_VALUES);
            coef = dF4_2dvti1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos12,&coef,INSERT_VALUES);
            coef = dF4_2dvti0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos13,&coef,INSERT_VALUES);
            coef = dF4_2dvtk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos14,&coef,INSERT_VALUES);
            coef = dF4_2dvtk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos15,&coef,INSERT_VALUES);
            coef = dF4_2dvr(i,j,k);
            MatSetValues(A,1,&pos,1,&pos16,&coef,INSERT_VALUES);
            coef = dF4_2dvri1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos17,&coef,INSERT_VALUES);
            coef = dF4_2dvri0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos18,&coef,INSERT_VALUES);
            coef = dF4_2dvp(i,j,k);
            MatSetValues(A,1,&pos,1,&pos19,&coef,INSERT_VALUES);
            coef = dF4_2dvpk1(i,j,k);
            MatSetValues(A,1,&pos,1,&pos20,&coef,INSERT_VALUES);
            coef = dF4_2dvpk0(i,j,k);
            MatSetValues(A,1,&pos,1,&pos21,&coef,INSERT_VALUES);
          }
        }
      }
    }
  }
 if(rank == 4){
     PetscInt I = npsi-1,J = n2th-1,K = nphi;
       for(i = 2;i < I; i++){
         for(j = 1;j <= J; j++){
           for(k = 1;k <= K; k++){
              PetscInt pos0=Istart2+(i-1)*J*K+(j-1)*K+k-1;//Br(i,j,k)
              PetscInt pos1=Istart3+(i-1)*J*K+(j-1)*K+k-1;//Bt(i,j,k)
              PetscInt pos2=Istart4+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
              PetscInt pos3=Istart5+(i-1)*J*K+(j-1)*K+k-1;//vr(i,j,k)
              PetscInt pos4=Istart6+(i-1)*J*K+(j-1)*K+k-1;//vt(i,j,k)
              PetscInt pos5=Istart7+(i-1)*J*K+(j-1)*K+k-1;//vp(i,j,k)
              PetscInt pos6=Istart2+i*J*K+(j-1)*K+k-1    ;//Br(i+1,j,k)
              PetscInt pos7=Istart4+i*J*K+(j-1)*K+k-1    ;//Bp(i+1,j,k)
              PetscInt pos8=Istart5+i*J*K+(j-1)*K+k-1    ;//vr(i+1,j,k)
              PetscInt pos9=Istart7+i*J*K+(j-1)*K+k-1    ;//vp(i+1,j,k)
              PetscInt pos10=Istart2+(i-2)*J*K+(j-1)*K+k-1;//Br(i-1,j,k)
              PetscInt pos11=Istart4+(i-2)*J*K+(j-1)*K+k-1;//Bp(i-1,j,k)
              PetscInt pos12=Istart5+(i-2)*J*K+(j-1)*K+k-1;//vr(i-1,j,k)
              PetscInt pos13=Istart7+(i-2)*J*K+(j-1)*K+k-1;//vp(i-1,j,k)
              PetscInt pos14=Istart3+(i-1)*J*K+(j-1)*K+k  ;//Bt(i,j,k+1)
              PetscInt pos15=Istart4+(i-1)*J*K+(j-1)*K+k  ;//Bp(i,j,k+1)
              PetscInt pos16=Istart6+(i-1)*J*K+(j-1)*K+k  ;//vt(i,j,k+1)
              PetscInt pos17=Istart7+(i-1)*J*K+(j-1)*K+k  ;//vp(i,j,k+1)
              PetscInt pos18=Istart3+(i-1)*J*K+(j-1)*K+k-2;//Bt(i,j,k-1)
              PetscInt pos19=Istart4+(i-1)*J*K+(j-1)*K+k-2;//Bp(i,j,k-1)
              PetscInt pos20=Istart6+(i-1)*J*K+(j-1)*K+k-2;//vt(i,j,k-1)
              PetscInt pos21=Istart7+(i-1)*J*K+(j-1)*K+k-2;//vp(i,j,k-1)
              if (k == K) {
                pos14=Istart3+(i-1)*J*K+(j-1)*K    ;//Bt(i,j,k+1)
                pos15=Istart4+(i-1)*J*K+(j-1)*K    ;//Bp(i,j,k+1)
                pos16=Istart6+(i-1)*J*K+(j-1)*K    ;//vt(i,j,k+1)
                pos17=Istart7+(i-1)*J*K+(j-1)*K    ;//vp(i,j,k+1)
              }
              if (k == 1) {
                pos18=Istart3+(i-1)*J*K+(j-1)*K+K-1;//Bt(i,j,k-1)
                pos19=Istart4+(i-1)*J*K+(j-1)*K+K-1;//Bp(i,j,k-1)
                pos20=Istart6+(i-1)*J*K+(j-1)*K+K-1;//vt(i,j,k-1)
                pos21=Istart7+(i-1)*J*K+(j-1)*K+K-1;//vp(i,j,k-1)
              }

             coef = dF4_3dBr(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBt(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBp(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvr(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos3,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvt(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos4,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvp(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos5,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBri1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos6,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBpi1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos7,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvri1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos8,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvpi1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos9,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBri0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos10,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBpi0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos11,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvri0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos12,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvpi0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos13,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBtk1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos14,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBpk1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos15,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvtk1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos16,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvpk1(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos17,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBtk0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos18,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dBpk0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos19,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvtk0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos20,&coef,INSERT_VALUES);CHKERRQ(ierr);
             coef = dF4_3dvpk0(i,j,k);
             ierr = MatSetValues(A,1,&pos0,1,&pos21,&coef,INSERT_VALUES);CHKERRQ(ierr);
           }
        }
      }
      i = 1;
      for(j=1;j<=J;j++){
        for(k=1;k<=K;k++){
          PetscInt pos0 = Istart4+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
          coef = 1.0;
          ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
          PetscInt pos1 = Istart4+i*J*K+(j-1)*K+k-1;//Bp(i+1,j,k)
          coef = -4./3.;
          ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
          PetscInt pos2 = Istart4+(i+1)*J*K+(j-1)*K+k-1;//Bp(i+2,j,k)
          coef = 1./3.;
          ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
     i = I;
     for(j=1;j<=J;j++){
       for(k=1;k<=K;k++){
         PetscInt pos0 = Istart4+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
         coef = 1.0;
         ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos1 = Istart4+(i-2)*J*K+(j-1)*K+k-1;//Bp(i-1,j,k)
         coef = -4./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos2 = Istart4+(i-3)*J*K+(j-1)*K+k-1;//Bp(i-2,j,k)
         coef = 1./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
       }
     }
  }
  if(rank == 5){
    PetscInt I = npsi-1,J = n2th-1,K = nphi;
     for(i=2;i<I;i++){
       for(j=1;j<=J;j++){
         for(k=1;k<=K;k++){
           PetscInt pos0 = Istart5+(i-1)*J*K+(j-1)*K+k-1;//vr(i,j,k)
           PetscInt pos1 = Istart6+(i-1)*J*K+(j-1)*K+k-1;//vt(i,j,k)
           PetscInt pos2 = Istart7+(i-1)*J*K+(j-1)*K+k-1;//vp(i,j,k)
           PetscInt pos3 = Istart0+(i-1)*J*K+(j-1)*K+k-1;//rho(i,j,k)
           PetscInt pos4 = Istart3+(i-1)*J*K+(j-1)*K+k-1;//Bt(i,j,k)
           PetscInt pos5 = Istart4+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
           PetscInt pos6 = Istart5+i*J*K+(j-1)*K+k-1    ;//vr(i+1,j,k)
           PetscInt pos7 = Istart5+(i-2)*J*K+(j-1)*K+k-1;//vr(i-1.j.k)
           PetscInt pos8 = Istart5+(i-1)*J*K+j*K+k-1    ;//vr(i,j+1,k)
           PetscInt pos9 = Istart5+(i-1)*J*K+(j-2)*K+k-1;//vr(i,j-1,k)
           PetscInt pos10= Istart5+(i-1)*J*K+(j-1)*K+k  ;//vr(i,j,k+1)
           PetscInt pos11= Istart5+(i-1)*J*K+(j-1)*K+k-2;//vr(i,j,k-1)
           PetscInt pos12= Istart1+i*J*K+(j-1)*K+k-1    ;//P(i+1,j,k)
           PetscInt pos13= Istart1+(i-2)*J*K+(j-1)*K+k-1;//P(i-1,j,k)
           PetscInt pos14= Istart2+(i-1)*J*K+j*K+k-1    ;//Br(i,j+1,k)
           PetscInt pos15= Istart2+(i-1)*J*K+(j-2)*K+k-1;//Br(i,j-1,k)
           PetscInt pos16= Istart2+(i-1)*J*K+(j-1)*K+k  ;//Br(i,j,k+1)
           PetscInt pos17= Istart2+(i-1)*J*K+(j-1)*K+k-2;//Br(i,j,k-1)
           PetscInt pos18= Istart3+i*J*K+(j-1)*K+k-1    ;//Bt(i+1,j,k)
           PetscInt pos19= Istart3+(i-2)*J*K+(j-1)*K+k-1;//Bt(i-1,j,k)
           PetscInt pos20= Istart4+i*J*K+(j-1)*K+k-1    ;//Bp(i+1,j,k)
           PetscInt pos21= Istart4+(i-2)*J*K+(j-1)*K+k-1;//Bp(i-1,j,k)

           if (j == 1){
            pos9 = Istart5+(i-1)*J*K+(J-1)*K+k-1;//vr(i,j-1,k)
            pos15= Istart2+(i-1)*J*K+(J-1)*K+k-1;//Br(i,j-1,k)
           }
           if (j == J){
            pos8 = Istart5+(i-1)*J*K+k-1    ;//vr(i,j+1,k)
            pos14= Istart2+(i-1)*J*K+k-1    ;//Br(i,j+1,k)
           }
           if (k == 1){
            pos11= Istart5+(i-1)*J*K+(j-1)*K+K-1;//vr(i,j,k-1)
            pos17= Istart2+(i-1)*J*K+(j-1)*K+K-1;//Br(i,j,k-1)
           }
           if (k == K){
            pos10= Istart5+(i-1)*J*K+(j-1)*K  ;//vr(i,j,k+1)
            pos16= Istart2+(i-1)*J*K+(j-1)*K  ;//Br(i,j,k+1)
           }


           coef = dF3_1dvr(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvt(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvp(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1drho(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos3,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBt(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos4,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBp(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos5,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvri1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos6,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvri0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos7,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvrj1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos8,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvrj0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos9,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvrk1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos10,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dvrk0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos11,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dPi1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos12,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dPi0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos13,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBrj1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos14,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBrj0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos15,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBrk1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos16,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBrk0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos17,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBti1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos18,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBti0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos19,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBpi1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos20,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_1dBpi0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos21,&coef,INSERT_VALUES);CHKERRQ(ierr);
         }
       }
     }
     i = 1;
     for(j=1;j<=J;j++){
       for(k=1;k<=K;k++){
         PetscInt pos0 = Istart5+(i-1)*J*K+(j-1)*K+k-1;//vr(i,j,k)
         coef = 1.0;
         ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos1 = Istart5+i*J*K+(j-1)*K+k-1;//vr(i+1,j,k)
         coef = -4./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos2 = Istart5+(i+1)*J*K+(j-1)*K+k-1;//vr(i+2,j,k)
         coef = 1./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
       }
     }
    i = I;
    for(j=1;j<=J;j++){
      for(k=1;k<=K;k++){
        PetscInt pos0 = Istart4+(i-1)*J*K+(j-1)*K+k-1;//vr(i,j,k)
        coef = 1.0;
        ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
        PetscInt pos1 = Istart4+(i-2)*J*K+(j-1)*K+k-1;//vr(i-1,j,k)
        coef = -4./3.;
        ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
        PetscInt pos2 = Istart4+(i-3)*J*K+(j-1)*K+k-1;//vr(i-2,j,k)
        coef = 1./3.;
        ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
      }
    }

  }
  if(rank == 6){
     PetscInt I = npsi-1,J = n2th-1,K = nphi;
     for(i=2;i<I;i++){
       for(j=1;j<=J;j++){
         for(k=1;k<=K;k++){
           PetscInt pos0 = Istart5+(i-1)*J*K+(j-1)*K+k-1;//vr(i,j,k)
           PetscInt pos1 = Istart6+(i-1)*J*K+(j-1)*K+k-1;//vt(i,j,k)
           PetscInt pos2 = Istart7+(i-1)*J*K+(j-1)*K+k-1;//vp(i,j,k)
           PetscInt pos3 = Istart0+(i-1)*J*K+(j-1)*K+k-1;//rho(i,j,k)
           PetscInt pos4 = Istart2+(i-1)*J*K+(j-1)*K+k-1;//Br(i,j,k)
           PetscInt pos5 = Istart3+(i-1)*J*K+(j-1)*K+k-1;//Bt(i,j,k)
           PetscInt pos6 = Istart4+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
           PetscInt pos7 = Istart6+i*J*K+(j-1)*K+k-1    ;//vt(i+1.j.k)
           PetscInt pos8 = Istart6+(i-2)*J*K+(j-1)*K+k-1;//vt(i-1,j,k)
           PetscInt pos9 = Istart6+(i-1)*J*K+j*K+k-1    ;//vt(i,j+1,k)
           PetscInt pos10= Istart6+(i-1)*J*K+(j-2)*K+k-1;//vt(i,j-1,k)
           PetscInt pos11= Istart6+(i-1)*J*K+(j-1)*K+k  ;//vt(i,j,k+1)
           PetscInt pos12= Istart6+(i-1)*J*K+(j-1)*K+k-2;//vt(i,j,k-1)
           PetscInt pos13= Istart1+(i-1)*J*K+j*K+k-1    ;//P(i,j+1,k)
           PetscInt pos14= Istart1+(i-1)*J*K+(j-2)*K+k-1;//P(i,j-1,k)
           PetscInt pos15= Istart2+(i-1)*J*K+j*K+k-1    ;//Br(i,j+1,k)
           PetscInt pos16= Istart2+(i-1)*J*K+(j-2)*K+k-1;//Br(i,j-1,k)
           PetscInt pos17= Istart3+i*J*K+(j-1)*K+k-1    ;//Bt(i+1,j,k)
           PetscInt pos18= Istart3+(i-2)*J*K+(j-1)*K+k-1;//Bt(i-1,j,k)
           PetscInt pos19= Istart3+(i-1)*J*K+(j-1)*K+k  ;//Bt(i,j,k+1)
           PetscInt pos20= Istart3+(i-1)*J*K+(j-1)*K+k-2;//Bt(i,j,k-1)
           PetscInt pos21= Istart4+(i-1)*J*K+j*K+k-1    ;//Bp(i,j+1,k)
           PetscInt pos22= Istart4+(i-1)*J*K+(j-2)*K+k-1;//Bp(i,j-1,k)

           if (j == 1){
             pos10= Istart6+(i-1)*J*K+(J-1)*K+k-1;//vt(i,j-1,k)
             pos14= Istart1+(i-1)*J*K+(J-1)*K+k-1;//P(i,j-1,k)
             pos16= Istart2+(i-1)*J*K+(J-1)*K+k-1;//Br(i,j-1,k)
             pos22= Istart4+(i-1)*J*K+(J-1)*K+k-1;//Bp(i,j-1,k)
           }
           if (j == J){
             pos9 = Istart6+(i-1)*J*K+k-1    ;//vt(i,j+1,k)
             pos13= Istart1+(i-1)*J*K+k-1    ;//P(i,j+1,k)
             pos15= Istart2+(i-1)*J*K+k-1    ;//Br(i,j+1,k)
             pos21= Istart4+(i-1)*J*K+k-1    ;//Bp(i,j+1,k)
           }
           if (k == 1){
             pos12= Istart6+(i-1)*J*K+(j-1)*K+K-1;//vt(i,j,k-1)
             pos20= Istart3+(i-1)*J*K+(j-1)*K+K-1;//Bt(i,j,k-1)
           }
           if (k ==K){
             pos11= Istart6+(i-1)*J*K+(j-1)*K  ;//vt(i,j,k+1)
             pos19= Istart3+(i-1)*J*K+(j-1)*K  ;//Bt(i,j,k+1)
           }

           coef = dF3_2dvr(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvt(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvp(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2drho(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos3,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBr(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos4,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBt(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos5,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBp(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos6,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvti1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos7,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvti0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos8,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvtj1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos9,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvtj0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos10,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvtk1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos11,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dvtk0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos12,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dPj1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos13,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dPj0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos14,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBrj1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos15,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBrj0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos16,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBti1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos17,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBti0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos18,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBtk1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos19,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBtk0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos20,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBpj1(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos21,&coef,INSERT_VALUES);CHKERRQ(ierr);
           coef = dF3_2dBpj0(i,j,k);
           ierr = MatSetValues(A,1,&pos0,1,&pos22,&coef,INSERT_VALUES);CHKERRQ(ierr);
         }
       }
     }
     i = 1;
     for(j=1;j<=J;j++){
       for(k=1;k<=K;k++){
         PetscInt pos0 = Istart6+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
         coef = 1.0;
         ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos1 = Istart6+i*J*K+(j-1)*K+k-1;//Bp(i+1,j,k)
         coef = -4./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos2 = Istart6+(i+1)*J*K+(j-1)*K+k-1;//Bp(i+2,j,k)
         coef = 1./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
       }
     }
    i = I;
    for(j=1;j<=J;j++){
      for(k=1;k<=K;k++){
        PetscInt pos0 = Istart6+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
        coef = 1.0;
        ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
        PetscInt pos1 = Istart6+(i-2)*J*K+(j-1)*K+k-1;//Bp(i-1,j,k)
        coef = -4./3.;
        ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
        PetscInt pos2 = Istart6+(i-3)*J*K+(j-1)*K+k-1;//Bp(i-2,j,k)
        coef = 1./3.;
        ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }
  if(rank == 7){
      PetscInt I = npsi-1,J = n2th-1,K = nphi;
      for(i=2;i<I;i++){
        for(j=1;j<=J;j++){
          for(k=1;k<=K;k++){
            PetscInt pos0 = Istart5+(i-1)*J*K+(j-1)*K+k-1;//vr(i,j,k)
            PetscInt pos1 = Istart6+(i-1)*J*K+(j-1)*K+k-1;//vt(i,j,k)
            PetscInt pos2 = Istart7+(i-1)*J*K+(j-1)*K+k-1;//vp(i,j,k)
            PetscInt pos3 = Istart0+(i-1)*J*K+(j-1)*K+k-1;//rho(i,j,k)
            PetscInt pos4 = Istart2+(i-1)*J*K+(j-1)*K+k-1;//Br(i,j,k)
            PetscInt pos5 = Istart3+(i-1)*J*K+(j-1)*K+k-1;//Bt(i,j,k)
            PetscInt pos6 = Istart4+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
            PetscInt pos7 = Istart7+i*J*K+(j-1)*K+k-1    ;//vp(i+1.j.k)
            PetscInt pos8 = Istart7+(i-2)*J*K+(j-1)*K+k-1;//vp(i-1,j,k)
            PetscInt pos9 = Istart7+(i-1)*J*K+j*K+k-1    ;//vp(i,j+1,k)
            PetscInt pos10= Istart7+(i-1)*J*K+(j-2)*K+k-1;//vp(i,j-1,k)
            PetscInt pos11= Istart7+(i-1)*J*K+(j-1)*K+k  ;//vp(i,j,k+1)
            PetscInt pos12= Istart7+(i-1)*J*K+(j-1)*K+k-2;//vp(i,j,k-1)
            PetscInt pos13= Istart1+(i-1)*J*K+(j-1)*K+k  ;//P(i,j,k+1)
            PetscInt pos14= Istart1+(i-1)*J*K+(j-1)*K+k-2;//P(i,j,k-1)
            PetscInt pos15= Istart2+(i-1)*J*K+(j-1)*K+k  ;//Br(i,j,k+1)
            PetscInt pos16= Istart2+(i-1)*J*K+(j-1)*K+k-2;//Br(i,j,k-1)
            PetscInt pos17= Istart3+(i-1)*J*K+(j-1)*K+k  ;//Bt(i,j,k+1)
            PetscInt pos18= Istart3+(i-1)*J*K+(j-1)*K+k-2;//Bt(i,j,k-1)
            PetscInt pos19= Istart4+i*J*K+(j-1)*K+k-1    ;//Bp(i+1,j,k)
            PetscInt pos20= Istart4+(i-2)*J*K+(j-1)*K+k-1;//Bp(i-1,j,k)
            PetscInt pos21= Istart4+(i-1)*J*K+j*K+k-1    ;//Bp(i,j+1,k)
            PetscInt pos22= Istart4+(i-1)*J*K+(j-2)*K+k-1;//Bp(i,j-1,k)

            if (j == 1){
              pos10= Istart7+(i-1)*J*K+(J-1)*K+k-1;//vp(i,j-1,k)
              pos22= Istart4+(i-1)*J*K+(J-1)*K+k-1;//Bp(i,j-1,k)
            }
            if (j == J){
              pos9 = Istart7+(i-1)*J*K+k-1    ;//vp(i,j+1,k)
              pos21= Istart4+(i-1)*J*K+k-1    ;//Bp(i,j+1,k)
            }
            if (k == 1){
              pos12= Istart7+(i-1)*J*K+(j-1)*K+K-1;//vp(i,j,k-1)
              pos14= Istart1+(i-1)*J*K+(j-1)*K+K-1;//P(i,j,k-1)
              pos16= Istart2+(i-1)*J*K+(j-1)*K+K-1;//Br(i,j,k-1)
              pos18= Istart3+(i-1)*J*K+(j-1)*K+K-1;//Bt(i,j,k-1)
            }
            if (k == K){
              pos11= Istart7+(i-1)*J*K+(j-1)*K  ;//vp(i,j,k+1)
              pos13= Istart1+(i-1)*J*K+(j-1)*K  ;//P(i,j,k+1)
              pos15= Istart2+(i-1)*J*K+(j-1)*K  ;//Br(i,j,k+1)
              pos17= Istart3+(i-1)*J*K+(j-1)*K  ;//Bt(i,j,k+1)
            }



            coef = dF3_3dvr(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvt(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvp(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3drho(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos3,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBr(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos4,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBt(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos5,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBp(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos6,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvpi1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos7,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvpi0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos8,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvpj1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos9,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvpj0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos10,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvpk1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos11,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dvpk0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos12,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dPk1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos13,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dPk0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos14,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBrk1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos15,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBrk0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos16,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBtk1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos17,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBtk0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos18,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBpi1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos19,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBpi0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos20,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBpj1(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos21,&coef,INSERT_VALUES);CHKERRQ(ierr);
            coef = dF3_3dBpj0(i,j,k);
            ierr = MatSetValues(A,1,&pos0,1,&pos22,&coef,INSERT_VALUES);CHKERRQ(ierr);
          }
        }
      }
      i = 1;
      for(j=1;j<=J;j++){
        for(k=1;k<=K;k++){
          PetscInt pos0 = Istart7+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
          coef = 1.0;
          ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
          PetscInt pos1 = Istart7+i*J*K+(j-1)*K+k-1;//Bp(i+1,j,k)
          coef = -4./3.;
          ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
          PetscInt pos2 = Istart7+(i+1)*J*K+(j-1)*K+k-1;//Bp(i+2,j,k)
          coef = 1./3.;
          ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
        }
      }
     i = I;
     for(j=1;j<=J;j++){
       for(k=1;k<=K;k++){
         PetscInt pos0 = Istart7+(i-1)*J*K+(j-1)*K+k-1;//Bp(i,j,k)
         coef = 1.0;
         ierr = MatSetValues(A,1,&pos0,1,&pos0,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos1 = Istart7+(i-2)*J*K+(j-1)*K+k-1;//Bp(i-1,j,k)
         coef = -4./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos1,&coef,INSERT_VALUES);CHKERRQ(ierr);
         PetscInt pos2 = Istart7+(i-3)*J*K+(j-1)*K+k-1;//Bp(i-2,j,k)
         coef = 1./3.;
         ierr = MatSetValues(A,1,&pos0,1,&pos2,&coef,INSERT_VALUES);CHKERRQ(ierr);
       }
     }
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);


  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr);

  //初始化b向量
  if(rank == 0){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F1(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }
  if(rank == 1){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F2(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }
  if(rank == 2){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F3_1(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }
  if(rank == 3){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F3_2(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }
  if(rank == 4){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F3_3(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }
  if(rank == 5){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F4_1(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }
  if(rank == 6){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F4_2(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }
  if(rank == 7){
    for(i = 1;i <= npsi-1;i++){
      for(j = 1;j <= n2th-1;j++){
        for(k = 1;k <= nphi;k++){
          PetscInt pos = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
          coef = F4_3(i,j,k);
          VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
        }
      }
    }
  }


  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-view_exact_sol",&flg,NULL);CHKERRQ(ierr);
  if (flg) {ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-2/(n+1),1.e-50,PETSC_DEFAULT,\
                          PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  for(n = 1;n <= nstep;i++){
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    KSPSolve(ksp,b,x);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    if(rank == 0){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            rho[i][j][k] = rho[i][j][k] + coef;
          }
        }
      }
    }
    if(rank == 1){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            P[i][j][k] = P[i][j][k] + coef;
          }
        }
      }
    }
    if(rank == 2){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            Br[i][j][k] = Br[i][j][k] + coef;
          }
        }
      }
    }
    if(rank == 3){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            Bt[i][j][k] = Bt[i][j][k] + coef;
          }
        }
      }
    }
    if(rank == 4){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            Bp[i][j][k] = Bp[i][j][k] + coef;
          }
        }
      }
    }
    if(rank == 5){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            vr[i][j][k] = vr[i][j][k] + coef;
          }
        }
      }
    }
    if(rank == 6){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            vt[i][j][k] = vt[i][j][k] + coef;
          }
        }
      }
    }
    if(rank == 7){
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            VecGetValues(x,1,&pos,&coef);
            vp[i][j][k] = vp[i][j][k] + coef;
          }
        }
      }
    }
    VecDestroy(&x);
    //重置b同时记录数据
    if(rank == 0){
      char fname1[] = "rho.dat";
      PetscFOpen(PETSC_COMM_SELF,fname1,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)rho[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart0 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F1(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }
    if(rank == 1){
      char fname2[] = "P.dat";
      PetscFOpen(PETSC_COMM_SELF,fname2,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)P[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart1 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F2(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }
    if(rank == 2){
      char fname3[] = "Br.dat";
      PetscFOpen(PETSC_COMM_SELF,fname3,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)Br[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart2 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F3_1(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }
    if(rank == 3){
      char fname4[] = "Bt.dat";
      PetscFOpen(PETSC_COMM_SELF,fname4,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)Bt[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart3 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F3_2(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }
    if(rank == 4){
      char fname5[] = "Bp.dat";
      PetscFOpen(PETSC_COMM_SELF,fname5,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)Bp[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart4 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F3_3(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }
    if(rank == 5){
      char fname6[] = "vr.dat";
      PetscFOpen(PETSC_COMM_SELF,fname6,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)vr[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart5 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F4_1(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }
    if(rank == 6){
      char fname7[] = "vt.dat";
      PetscFOpen(PETSC_COMM_SELF,fname7,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)vt[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart6 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F4_2(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }
    if(rank == 7){
      char fname8[] = "vp.dat";
      PetscFOpen(PETSC_COMM_SELF,fname8,"w",&fp);
      if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
      for(i = 1;i <= npsi-1;i++){
        for(j = 1;j <= n2th-1;j++){
          PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)vp[i][j][1]);
          for(k = 1;k <= nphi;k++){
            PetscInt pos = Istart7 + (i-1)*(n2th-1)*(nphi) + (j-1)*(nphi) + k-1;
            coef = F4_3(i,j,k);
            VecSetValues(b,1,&pos,&coef,INSERT_VALUES);
          }
        }
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
      PetscFClose(PETSC_COMM_SELF,fp);
    }

  }

  ierr = PetscFinalize( );
  return 0;
}
