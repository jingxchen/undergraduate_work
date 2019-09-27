/*
f=rho,P,v,B;t=theta;p=phi;
fri0=fri-1jk;fr=frijk;fri1=fri+1jk;
frj0=frij-1k;fr=frijk;frj1=frij+1k;
frk0=frijk-1;fr=frijk;frk1=frijk+1;
*/



--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

PetscScalar dF1drho(PetscScalar dt,PetscScalar r ,PetscScalar a,PetscScalar theta,PetscScalar vr,PetscScalar dr,PetscScalar vri1,PetscScalar vri0,PetscScalar vt,
                    PetscScalar vtj1,PetscScalar vtj0,PetscScalar dtheta,PetscScalar vpk1,PetscScalar vpk0,PetscScalar dphi)
{
       PetscScalar s;
       s=1/dt+(a+2*r*cos(theta))/(r*(a+r*cos(theta)))*vr+(vri1-vri0)/(2*dr)-sin(theta)/(a+r*cos(theta))*vt+1/r*(vtj1-vtj0)/(2*dtheta)+(vpk1-vpk0)/(2*dphi*(a+cos(theta)));
       return s;
}

PetscScalar dF1drhoi0(PetscScalar vr,PetscScalar dr){
  PetscScalar s;
  s=-vr/(2*dr);
  return s;
}

PetscScalar dF1drhoi1(PetscScalar vr ,PetscScalar dr)
{
  PetscScalar s;
  s=vr/(2*dr);
  return s;
}

PetscScalar dF1drhoj0(PetscScalar vt,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s=-vt/(2*r*dtheta);
  return s;
}

PetscScalar dF1drhoj1(PetscScalar vt,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s=vt/(2*r*dtheta);
  return s;
}

PetscScalar dF1drhok0(PetscScalar vp,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = -vp/((a+r*cos(theta))*dphi*2);
  return s;
}

PetscScalar dF1drhok1(PetscScalar vp,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = vp/((a+r*cos(theta))*dphi*2);
  return s;
}

PetscScalar dF1dvri0(PetscScalar rho ,PetscScalar dr)
{
  PetscScalar s;
  s= -rho/(2*dr);
  return s;
}

PetscScalar dF1dvri1(PetscScalar rho ,PetscScalar dr)
{
  PetscScalar s;
  s= rho/(2*dr);
  return s;
}

PetscScalar dF1dvtj0(PetscScalar rho, PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
   s = -rho/(2*r*dtheta);
   return s;
}

PetscScalar dF1dvtj1(PetscScalar rho, PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
   s = rho/(2*r*dtheta);
   return s;
}

PetscScalar dF1dvpk0(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s= -rho/((a+r*cos(theta))*2*dphi);
  return s;
}

PetscScalar dF1dvpk1(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s= rho/((a+r*cos(theta))*2*dphi);
  return s;
}

PetscScalar dF1dvr(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar rho,PetscScalar rhoi1,PetscScalar rhoi0,PetscScalar dr)
{
  PetscScalar s;
  s= (a+r*cos(theta)*2)/(r*(a+r*cos(theta)))*rho+(rhoi1-rhoi0)/(2*dr);
  return s;
}

PetscScalar dF1dvt(PetscScalar theta,PetscScalar a,PetscScalar r,PetscScalar rho,PetscScalar rhoj1,PetscScalar rhoj0,PetscScalar dtheta)
{
  PetscScalar s;
  s = -sin(theta)/(a+r*cos(theta))*rho+1/r*(rhoj1-rhoj0)/(2*dtheta);
  return s;
}

PetscScalar dF1dvp(PetscScalar rhok1,PetscScalar rhok0,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = (rhok1-rhok0)/((a+r*cos(theta))*2*dphi);
  return s;
}
-----------------------------------------------------------------------------------------------------------------------------------------







PetscScalar dF2dP(PetscScalar dt,PetscScalar r,PetscScalar a,PetscScalar theta,PetscScalar vri1,PetscScalar vri0,PetscScalar dr,PetscScalar vr,
                  PetscScalar vtj1,PetscScalar vtj0,PetscScalar dtheta,PetscScalar vt,PetscScalar vpk1,PetscScalar vpk0,PetscScalar dphi,PetscScalar gamma)
{
	PetscScalar s;
	s=-1/dt-gamma/(r*(a+r*cos(theta)))*(r*(a+r*cos(theta))*(vri1-vri0)/(2*dr)+(a+2*r*cos(theta))*vr+(a+r*cos(theta))*(vtj1-vtj0)/(2*dtheta)-r*sin(theta)*vt+r*(vpk1-vpk0)/(2*dphi));
	return s;
}

PetscScalar dF2dvr(PetscScalar dr,PetscScalar Pi1,PetscScalar Pi0,PetscScalar gamma,PetscScalar r,PetscScalar a,PetscScalar theta,PetscScalar P)
{
	PetscScalar s;
	s = -(Pi1-Pi0)/(2*dr)-gamma*P/(r*(a+r*cos(theta)))*(a+2*r*cos(theta));
	return s;
}

PetscScalar dF2dvt(PetscScalar r,PetscScalar Pj1,PetscScalar Pj0,PetscScalar dtheta,PetscScalar gamma,PetscScalar P,PetscScalar theta,PetscScalar a)
{
	PetscScalar s;
	s=-1/r*(Pj1-Pj0)/(2*dtheta)+gamma*P/(a+r*cos(theta))*sin(theta);
	return s;
}

PetscScalar dF2dvp(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Pk1,PetscScalar Pk0,PetscScalar dphi)
{
	PetscScalar s;
	s=-1/(a+r*cos(theta))*(Pk1-Pk0)/(2*dphi);
	return s;
}

PetscScalar dF2dPi1(PetscScalar vr,PetscScalar dr)
{
	PetscScalar s;
	s=-vr/(2*dr);
	return s;
}

PetscScalar dF2dPi0(PetscScalar vr,PetscScalar dr)
{
	PetscScalar s;
	s=vr/(2*dr);
	return s;
}

PetscScalar dF2dPj1(PetscScalar r,PetscScalar vt,PetscScalar dtheta)
{
	PetscScalar s;
	s=-vt/(2*r*dtheta);
	return s;
}

PetscScalar dF2dPj0(PetscScalar r,PetscScalar vt,PetscScalar dtheta)
{
	PetscScalar s;
	s=vt/(2*r*dtheta);
	return s;
}

PetscScalar dF2dPk1(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
	PetscScalar s;
	s=-1/(a+r*cos(theta))*vp/(2*dphi);
	return s;
}

PetscScalar dF2dPk0(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
	PetscScalar s;
	s=1/(a+r*cos(theta))*vp/(2*dphi);
	return s;
}

PetscScalar dF2dvri1(PetscScalar gamma,PetscScalar P,PetscScalar dr)
{
	PetscScalar s;
	s=-gamma*P/(2*dr);
	return s;
}

PetscScalar dF2dvri0(PetscScalar gamma,PetscScalar P,PetscScalar dr)
{
	PetscScalar s;
	s=gamma*P/(2*dr);
	return s;
}

PetscScalar dF2dvtj1(PetscScalar gamma,PetscScalar P,PetscScalar r,PetscScalar dtheta)
{
	PetscScalar s;
	s=-gamma*P/(r*2*dtheta);
	return s;
}

PetscScalar dF2dvtj0(PetscScalar gamma,PetscScalar P,PetscScalar r,PetscScalar dtheta)
{
	PetscScalar s;
	s=gamma*P/(r*2*dtheta);
	return s;
}

PetscScalar dF2dvpk1(PetscScalar gamma,PetscScalar P,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
	PetscScalar s;
	s=-gamma*P/((a+r*cos(theta))*2*dphi);
	return s;
}

PetscScalar dF2dvpk0(PetscScalar gamma,PetscScalar P,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
	PetscScalar s;
	s=gamma*P/((a+r*cos(theta))*2*dphi);
	return s;
}
------------------------------------------------------------------------------------------------------------------------------------------------------------------









PetscScalar dF3_1dvr(PetscScalar dt,PetscScalar dr,PetscScalar vri1,PetscScalar vri0)
{
  PetscScalar s;
  s = -1/dt-(vri1-vri0)/(2*dr);
  return s;
}

PetscScalar dF3_1dvt(PetscScalar r,PetscScalar vrj1,PetscScalar vrj0,PetscScalar dtheta)
{
  PetscScalar s;
  s = -1/r*(vrj1-vrj0)/(2*dtheta);
  return s;
}

PetscScalar dF3_1dvp(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,PetscScalar vrk1,PetscScalar vrk0)
{
  PetscScalar s;
  s = -1/(a+r*cos(theta))*(vrk1-vrk0)/(2*dphi);
  return s;
}

PetscScalar dF3_1drho(PetscScalar rho,PetscScalar dr,PetscScalar Pri1,PetscScalar Pri0,PetscScalar Bp,PetscScalar a,PetscScalar r,PetscScalar theta,
                      PetscScalar Brk1,PetscScalar Brk0,PetscScalar dphi,PetscScalar Bpi1,PetscScalar Bpi0,PetscScalar Bt,PetscScalar Bti1,PetscScalar Bti0,
                      PetscScalar Brj1,PetscScalar Brj0,PetscScalar dtheta)
{
  PetscScalar s;
  s = -1/(rho*rho)*(-(Pri1-Pri0)/(2*dr)+Bp/(a+r*cos(theta))*(Brk1-Brk0)/(2*dphi)-Bp*Bp*cos(theta)/(a+r*cos(theta))-Bp*(Bpi1-Bpi0)/(2*dr)-Bt*Bt/r-Bt*(Bti1-Bti0)/(2*dr)+(Brj1-Brj0)*Bt/(2*r*dtheta));
  return s;
}

PetscScalar dF3_1dBt(PetscScalar Bt,PetscScalar rho,PetscScalar r,PetscScalar Bti1,PetscScalar Bti0,PetscScalar dr,PetscScalar Brj1,PetscScalar Brj0,PetscScalar dtheta)
{
  PetscScalar s;
  s = -2*Bt/(rho*r)-(Bti1-Bti0)/(rho*2*dr)+(Brj1-Brj0)/(2*r*dtheta*rho);
  return s;
}

PetscScalar dF3_1dBp(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,PetscScalar Brk1,PetscScalar Brk0,PetscScalar dr,PetscScalar Bpi1,PetscScalar Bpi0,
                     PetscScalar Bp)
{
  PetscScalar s;
  s = (Brk1-Brk0)/(rho*(a+r*cos(theta))*2*dphi)-2*Bp*cos(theta)/((a+r*cos(theta))*rho-(Bpi1-Bpi0)/(rho*2*dr);
  return s;
}

PetscScalar dF3_1dvri1(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = -vr/(2*dr);
  return s;
}

PetscScalar dF3_1dvri0(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = vr/(2*dr);
  return s;
}

PetscScalar dF3_1dvrj1(PetscScalar vt,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = -vt/(r*2*dtheta);
  return s;
}

PetscScalar dF3_1dvrj0(PetscScalar vt,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = vt/(r*2*dtheta);
  return s;
}

PetscScalar dF3_1dvrk1(PetscScalar vp,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = -vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_1dvrk0(PetscScalar vp,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = vp/(2*dphi*(a+r*cos(theta)));
  return s;
}


PetscScalar dF3_1dPr1(PetscScalar rho,PetscScalar dr)
{
  PetscScalar s;
  s = -1/(rho*2*dr);
  return s;
}

PetscScalar dF3_1dPr0(PetscScalar rho,PetscScalar dr)
{
  PetscScalar s;
  s = 1/(rho*2*dr);
  return s;
}

PetscScalar dF3_1dBrj1(PetscScalar rho,PetscScalar r,PetscScalar dtheta,PetscScalar Bt)
{
  PetscScalar s;
  s = Bt/(rho*r*2*dtheta);
  return s;
}

PetscScalar dF3_1dBrj0(PetscScalar rho,PetscScalar r,PetscScalar dtheta,PetscScalar Bt)
{
  PetscScalar s;
  s = -Bt/(rho*r*2*dtheta);
  return s;
}

PetscScalar dF3_1dBrk1(PetscScalar Bp,PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = Bp/(rho*(a+r*cos(theta)*2*dphi));
  return s;
}

PetscScalar dF3_1dBrk0(PetscScalar Bp,PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = -Bp/(rho*(a+r*cos(theta)*2*dphi));
  return s;
}

PetscScalar dF3_1dBti1(PetscScalar Bt,PetscScalar dr,PetscScalar rho)
{
  PetscScalar s;
  s = -Bt/(2*dr*rho);
  return s;
}

PetscScalar dF3_1dBti0(PetscScalar Bt,PetscScalar dr,PetscScalar rho)
{
  PetscScalar s;
  s = Bt/(2*dr*rho);
  return s;
}


PetscScalar dF3_1dBpi1(PetscScalar Bp,PetscScalar dr,PetscScalar rho)
{
  PetscScalar s;
  s = -Bp/(2*dr*rho);
  return s;
}

PetscScalar dF3_1dBpi0(PetscScalar Bp,PetscScalar dr,PetscScalar rho)
{
  PetscScalar s;
  s = Bp/(2*dr*rho);
  return s;
}
-------------------------------------------------------------------------------------------------------------------------------------------------------------------











PetscScalar dF3_2dvr(PetscScalar vti1,PetscScalar vti0,PetscScalar dr)
{
  PetscScalar s;
  s = -(vti1-vti0)/(2*dr);
  return s;
}

PetscScalar dF3_2dvt(PetscScalar dt,PetscScalar r,PetscScalar dtheta,PetscScalar vtj1,PetscScalar vtj0)
{
  PetscScalar s;
  s = -1/dt-(vtj1-vtj0)/(r*2*dtheta);
  return s;
}

PetscScalar dF3_2dvp(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,PetscScalar vtk1,PetscScalar vtk0)
{
  PetscScalar s;
  s = -1/(a+r*cos(theta))*(vtk1-vtk0)/(2*dphi);
  return s;
}

PetscScalar dF3_2drho(PetscScalar rho,PetscScalar r,PetscScalar Ptj1,PetscScalar Ptj0,PetscScalar dtheta,PetscScalar Bt,PetscScalar Br,PetscScalar Bti1,PetscScalar Bti0,PetscScalar dr,PetscScalar Brj1,
                      PetscScalar Brj0,PetscScalar Bp,PetscScalar Bpj1,PetscScalar Bpj0,PetscScalar Btk1,PetscScalar Btk0)
{
  PetscScalar s;
  s = -1/(rho*rho)*(-(Ptj1-Ptj0)/(2*r*dtheta)+Bt*Br/r+Br*(Bti1-Bti0)/(2*dr)-Br/r*(Brj1-Brj0)/(2*dtheta)+(Bp*Bp)*sin(theta)/(a+r*cos(theta))-Bp/r*(Bpj1-Bpj0)/(2*dtheta)
      +Bp/(a+r*cos(theta))*(Btk1-Btk0)/(2*dphi));
  return s;
}

PetscScalar dF3_2dBr(PetscScalar rho,PetscScalar Bt,PetscScalar r,PetscScalar Bti1,PetscScalar Bti0,PetscScalar dr,PetscScalar Brj1,PetscScalar Brj0,PetscScalar dtheta)
{
  PetscScalar s;
  s = 1/rho*(Bt/r+(Bti1-Bti0)/(2*dr)-(Brj1-Brj0)/(r*2*dtheta));
  return s;
}

PetscScalar dF3_2dBt(PetscScalar rho,PetscScalar Br,PetscScalar r)
{
  PetscScalar s;
  s = Br/(rho*r);
  return s;
}

PetscScalar dF3_2dBp(PetscScalar rho,PetscScalar theta,PetscScalar a,PetscScalar r,PetscScalar Bp,PetscScalar Bpj1,PetscScalar Bpj0,PetscScalar dtheta,PetscScalar Btk1,PetscScalar Btk0,PetscScalar dphi)
{
  PetscScalar s;
  s = 1/rho*(2*sin(theta)*Bp/(a+r*cos(theta))-(Bpj1-Bpj0)/(2*r*dtheta)+(Btk1-Btk0)/(2*dphi*(a+r*cos(theta))));
  return s;
}

PetscScalar dF3_2dvti1(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = -vr/2*(dr);
  return s;
}

PetscScalar dF3_2dvti0(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = vr/2*(dr);
  return s;
}

PetscScalar dF3_2dvtj1(PetscScalar r,PetscScalar vt,PetscScalar dtheta)
{
  PetscScalar s;
  s = -vt/(r*2*dtheta);
  return s;
}

PetscScalar dF3_2dvtj0(PetscScalar r,PetscScalar vt,PetscScalar dtheta)
{
  PetscScalar s;
  s = vt/(r*2*dtheta);
  return s;
}

PetscScalar dF3_2dvtk1(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = -vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_2dvtk0(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_2dPtj1(PetscScalar rho,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = -1/(rho*r*2*dtheta);
  return s;
}

PetscScalar dF3_2dPtj0(PetscScalar rho,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = 1/(rho*r*2*dtheta);
  return s;
}

PetscScalar dF3_2dBrj1(PetscScalar rho,PetscScalar Br,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = -Br/(rho*r*2*dtheta);
  return s;
}

PetscScalar dF3_2dBrj0(PetscScalar rho,PetscScalar Br,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = Br/(rho*r*2*dtheta);
  return s;
}

PetscScalar dF3_2dBti1(PetscScalar rho,PetscScalar Br,PetscScalar dr)
{
  PetscScalar s;
  s = Br/(rho*2*dr);
  return s;
}

PetscScalar dF3_2dBti0(PetscScalar rho,PetscScalar Br,PetscScalar dr)
{
  PetscScalar s;
  s = -Br/(rho*2*dr);
  return s;
}


PetscScalar dF3_2dBtk1(PetscScalar Bp,PetscScalar dphi,PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = Bp/(2*dphi*rho*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_2dBtk0(PetscScalar Bp,PetscScalar dphi,PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = -Bp/(2*dphi*rho*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_2dBpj1(PetscScalar rho,PetscScalar Bp,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = -Bp/(rho*r*2*dtheta);
  return s;
}

PetscScalar dF3_2dBpj0(PetscScalar rho,PetscScalar Bp,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = Bp/(rho*r*2*dtheta);
  return s;
}
------------------------------------------------------------------------------------------------------------------------------------------------------------------










PetscScalar dF3_3dvr(PetscScalar dr,PetscScalar vpi1,PetscScalar vpi0)
{
  PetscScalar s;
  s = -(vpi1-vpi0)/(2*dr);
  return s;
}

PetscScalar dF3_3dvt(PetscScalar dtheta,PetscScalar r,PetscScalar vpj1,PetscScalar vpj0)
{
  PetscScalar s;
  s = -(vpj1-vpj0)/(r*2*dtheta);
  return s;
}

PetscScalar dF3_3dvp(PetscScalar dt,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,petsc vpk1,PetscScalar vpk0)
{
  PetscScalar s;
  s = -1/dt-1/(a+r*cos(theta))*(vpk1-vpk0)/(2*dphi);
  return s;
}

PetscScalar dF3_3drho(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,PetscScalar Ppk1,PetscScalar Ppk0,PetscScalar Bp,
                      PetscScalar Bt,PetscScalar dtheta,PetscScalar Bpj1,PetscScalar Bpj0,PetscScalar Btk1,PetscScalar Btk0,PetscScalar Br,PetscScalar Brk1,
                      PetscScalar Brk0,PetscScalar Bpi1,PetscScalar Bpi0,PetscScalar dr)
{
  PetscScalar s;
  s = -1/(rho*rho)*(-(Ppk1-Ppk0)/(2*dphi*(a+r*cos(theta)))-Bp*Bt*sin(theta)/(a+r*cos(theta))+Bt/r*(Bpj1-Bpj0)/(2*dtheta)-Bt/(a+r*cos(theta))*(Btk1-Btk0)/(2*dphi)
                    -Br/(a+r*cos(theta))*(Brk1-Brk0)/(2*dphi)+cos(theta)*Bp*Br/(a+r*cos(theta))+Br*(Bpi1-Bpi0)/(2*dr));
  return s;
}

PetscScalar dF3_3dBr(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,PetscScalar Brk1,PetscScalar Brk0,PetscScalar Bp,PetscScalar dr,
                     PetscScalar Bpi1,PetscScalar Bpi0)
{
  PetscScalar s;
  s = 1/rho*(-(Brk1-Brk0)/(2*dphi*(a+r*cos(theta)))+cos(theta)*Bp/(a+r*cos(theta))+(Bpi1-Bpi0)/(2*dr));
  return s;
}

PetscScalar dF3_3dBt(PetscScalar rho,PetscScalar Bp,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dtheta,PetscScalar Bpj1,PetscScalar Bpj0,PetscScalar dphi,
                     PetscScalar Btk1,PetscScalar Btk0)
{
  PetscScalar s;
  s = 1/rho*(-Bp*sin(theta)/(a+r*cos(theta))+(Bpj1-Bpj0)/(2*dtheta*r)-(Btk1-Btk0)/(2*dphi*(a+r*cos(theta))));
  return s;
}

PetscScalar dF3_3dBp(PetscScalar rho,PetscScalar Bt,PetscScalar theta,PetscScalar a,PetscScalar r,PetscScalar Br)
{
  PetscScalar s;
  s = 1/rho*(-Bt*sin(theta)/(a+r*cos(theta))+Br*cos(theta)/(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_3dvpi1(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = -vr/(2*dr);
  return s;
}

PetscScalar dF3_3dvpi0(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = vr/(2*dr);
  return s;
}

PetscScalar dF3_3dvpj1(PetscScalar vt,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = -vt/(r*2*dtheta);
  return s;
}

PetscScalar dF3_3dvpj0(PetscScalar vt,PetscScalar r,PetscScalar dtheta)
{
  PetscScalar s;
  s = vt/(r*2*dtheta);
  return s;
}

PetscScalar dF3_3dvpk1(PetscScalar vp,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = -vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_3dvpk0(PetscScalar vp,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_3dPpk1(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = -1/(rho*2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF3_3dPpk0(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi)
{
  PetscScalar s;
  s = 1/(rho*2*dphi*(a+r*cos(theta)));
  return s;
}


PetscScalar dF3_3dBrk1(PetscScalar rho,PetscScalar Br,PetscScalar dphi,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = -Br/((a+r*cos(theta))*rho*2*dphi);
  return s;
}

PetscScalar dF3_3dBrk0(PetscScalar rho,PetscScalar Br,PetscScalar dphi,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = Br/((a+r*cos(theta))*rho*2*dphi);
  return s;
}

PetscScalar dF3_3dBtk1(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Bt,PetscScalar dphi)
{
  PetscScalar s;
  s = -Bt/((a+r*cos(theta))*rho*2*dphi);
  return s;
}

PetscScalar dF3_3dBtk0(PetscScalar rho,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Bt,PetscScalar dphi)
{
  PetscScalar s;
  s = Bt/((a+r*cos(theta))*rho*2*dphi);
  return s;
}

PetscScalar dF3_3dBpi1(PetscScalar Br,PetscScalar dr,PetscScalar rho)
{
  PetscScalar s;
  s = Br/(2*dr*rho);
  return s;
}

PetscScalar dF3_3dBpi0(PetscScalar Br,PetscScalar dr,PetscScalar rho)
{
  PetscScalar s;
  s = -Br/(2*dr*rho);
  return s;
}

PetscScalar dF3_3dBpj1(PetscScalar Bt,PetscScalar r,PetscScalar dtheta,PetscScalar rho)
{
  PetscScalar s;
  s = Bt/(r*rho*2*dtheta);
  return s;
}
PetscScalar dF3_3dBpj0(PetscScalar Bt,PetscScalar r,PetscScalar dtheta,PetscScalar rho)
{
  PetscScalar s;
  s = -Bt/(r*rho*2*dtheta);
  return s;
}

------------------------------------------------------------------------------------------------------------------------------------------------------------------















PetscScalar dF4_1dBr(PetscScalar dt,PetscScalar r,PetscScalar dtheta,PetscScalar vtj1,PetscScalar vtj0,PetscScalar a,PetscScalar vpk1,PetscScalar vpk0,PetscScalar dphi)
{
  PetscScalar s;
  s = -1/dt-1/r*(vtj1-vtj0)/(2*dtheta)+sin(theta)/(a+r*cos(theta))*vt-(vpk1-vpk0)/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dBt(PetscScalar r,PetscScalar vrj1,PetscScalar vrj0,PetscScalar dtheta,PetscScalar theta ,PetscScalar a,PetscScalar vr)
{
  PetscScalar s;
  s = (vrj1-vrj0)/(2*r*dtheta)-sin(theta)*vr/(a+r*cos(theta));
  return s;
}

PetscScalar dF4_1dBp(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vrk1,PetscScalar vrk0,PetscScalar dphi)
{
  PetscScalar s;
  s = (vrk1-vrk0)/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dvr(PetscScalar r,PetscScalar Btj1,PetscScalar Btj0,PetscScalar dtheta ,PetscScalar theta,PetscScalar a,PetscScalar Bt,PetscScalar Bpk1,PetscScalar Bpk0,PetscScalar dphi)
{
  PetscScalar s;
  s = (Btj1-Btj0)/(2*r*dtheta)-sin(theta)*Bt/(a+r*cos(theta))+(Bpk1-Bpk0)/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dvt(PetscScalar r,PetscScalar dtheta,PetscScalar Brj1,PetscScalar Brj0,PetscScalar theta,PetscScalar a,PetscScalar Br)
{
  PetscScalar s;
  s = -(Brj1-Brj0)/(2*r*dtheta)+sin(theta)*Br/(a+r*cos(theta))
  return s;
}

PetscScalar dF4_1dvp(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Brk1,PetscScalar Brk0,PetscScalar dphi)
{
  PetscScalar s;
  s = -(Brk1-Brk0)/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dBrj1(PetscScalar r ,PetscScalar dtheta,PetscScalar vt)
{
  PetscScalar s;
  s = -vt/(2*r*dtheta);
  return s;
}

PetscScalar dF4_1dBrj0(PetscScalar r ,PetscScalar dtheta,PetscScalar vt)
{
  PetscScalar s;
  s = vt/(2*r*dtheta);
  return s;
}

PetscScalar dF4_1dBrk1(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = -vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dBrk0(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dBtj1(PetscScalar r,PetscScalar dtheta,PetscScalar vr)
{
  PetscScalar s;
  s = vr/(2*dtheta*r);
  return s;
}

PetscScalar dF4_1dBtj0(PetscScalar r,PetscScalar dtheta,PetscScalar vr)
{
  PetscScalar s;
  s = -vr/(2*dtheta*r);
  return s;
}

PetscScalar dF4_1dBpk1(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,PetscScalar vr)
{
  PetscScalar s;
  s = vr/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dBpk0(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar dphi,PetscScalar vr)
{
  PetscScalar s;
  s = -vr/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dvrj1(PetscScalar r,PetscScalar dtheta ,PetscScalar Bt)
{
  PetscScalar s;
  s = Bt/(r*2*dtheta);
  return s;
}

PetscScalar dF4_1dvrj0(PetscScalar r,PetscScalar dtheta ,PetscScalar Bt)
{
  PetscScalar s;
  s = -Bt/(r*2*dtheta);
  return s;
}

PetscScalar dF4_1dvrk1(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Bp,PetscScalar dphi)
{
  PetscScalar s;
  s = Bp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dvrk0(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Bp,PetscScalar dphi)
{
  PetscScalar s;
  s = -Bp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dvtj1(PetscScalar r,PetscScalar Br,PetscScalar dtheta)
{
  PetscScalar s;
  s = -Br/(2*dtheta*r);
  return s;
}

PetscScalar dF4_1dvtj0(PetscScalar r,PetscScalar Br,PetscScalar dtheta)
{
  PetscScalar s;
  s = Br/(2*dtheta*r);
  return s;
}

PetscScalar dF4_1dvpk1(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Br,PetscScalar dphi)
{
  PetscScalar s;
  s = -Br/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_1dvpk0(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Br,PetscScalar dphi)
{
  PetscScalar s;
  s = Br/(2*dphi*(a+r*cos(theta)));
  return s;
}
-----------------------------------------------------------------------------------------------------------------------------------------------------













PetscScalar dF4_2dBr(PetscScalar vti1,PetscScalar vti0,PetscScalar dr,PetscScalar theta,PetscScalar a,PetscScalar vt)
{
  PetscScalar s;
  s = (vti1-vti0)/(2*dr)+cos(theta)*vt/(a+r*cos(theta));
  return s;
}

PetscScalar dF4_2dBt(PetscScalar dt,PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vpk1,PetscScalar vpk0,PetscScalar dphi,PetscScalar vri1,PetscScalar vri0,PetscScalar dr,PetscScalar vr)
{
  PetscScalar s;
  s = -1/dt-(vpk1-vpk0)/(2*dphi*(a+r*cos(theta)))-(vri1-vri0)/(2*dr)-cos(theta)*dr/(a+r*cos(theta));
  return s;
}

PetscScalar dF4_2dBp(PetscScalar a,PetscScalar r,PetscScalar theta ,PetscScalar dphi,PetscScalar vtk1,PetscScalar vtk0)
{
  PetscScalar s;
  s = (vtk1-vtk0)/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dvr(PetscScalar Bti1,PetscScalar Bti0,PetscScalar dr,PetscScalar theta,PetscScalar a,PetscScalar r,PetscScalar Bt)
{
  PetscScalar s;
  s = -(Bti1-Bti0)/(2*dr)-cos(theta)*Bt/(a+r*cos(theta));
  return s;
}

PetscScalar dF4_2dvt(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Bpk1,PetscScalar Bpk0,PetscScalar dphi,PetscScalar Bri1,PetscScalar Bri0,PetscScalar dr,PetscScalar Br)
{
  PetscScalar s;
  s = (Bpk1-Bpk0)/(2*dphi*(a+r*cos(theta)))+(Bri1-Bri0)/(2*dr)+cos(theta)*Br/(a+r*cos(theta));
  return s;
}

PetscScalar dF4_2dvp(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar Btk1,PetscScalar Btk0,PetscScalar dphi)
{
  PetscScalar s;
  s = -(Btk1-Btk0)/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dBri1(PetscScalar vt,PetscScalar dr)
{
  PetscScalar s;
  s = vt/(2*dr);
  return s;
}

PetscScalar dF4_2dBri0(PetscScalar vt,PetscScalar dr)
{
  PetscScalar s;
  s = -vt/(2*dr);
  return s;
}

PetscScalar dF4_2dBti1(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = -vr/(2*dr);
  return s;
}

PetscScalar dF4_2dBti0(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = vr/(2*dr);
  return s;
}

PetscScalar dF4_2dBtk1(PetscScalar a,PetscScalar r ,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = -vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dBtk0(PetscScalar a,PetscScalar r ,PetscScalar theta,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = vp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dBpk1(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vt,PetscScalar dphi)
{
  PetscScalar s;
  s = vt/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dBpk0(PetscScalar a,PetscScalar r,PetscScalar theta,PetscScalar vt,PetscScalar dphi)
{
  PetscScalar s;
  s = -vt/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dvri1(PetscScalar Bt,PetscScalar dr)
{
  PetscScalar s;
  s = -Bt/(2*dr);
  return s;
}

PetscScalar dF4_2dvri0(PetscScalar Bt,PetscScalar dr)
{
  PetscScalar s;
  s = Bt/(2*dr);
  return s;
}

PetscScalar dF4_2dvti1(PetscScalar Br,PetscScalar dr)
{
  PetscScalar s;
  s = Br/(2*dr);
  return s;
}

PetscScalar dF4_2dvti0(PetscScalar Br,PetscScalar dr)
{
  PetscScalar s;
  s = -Br/(2*dr);
  return s;
}

PetscScalar dF4_2dvtk1(PetscScalar Bp,PetscScalar dphi,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = Bp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dvtk0(PetscScalar Bp,PetscScalar dphi,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = -Bp/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dvpk1(PetscScalar Bt,PetscScalar dphi,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = -Bt/(2*dphi*(a+r*cos(theta)));
  return s;
}

PetscScalar dF4_2dvpk0(PetscScalar Bt,PetscScalar dphi,PetscScalar a,PetscScalar r,PetscScalar theta)
{
  PetscScalar s;
  s = Bt/(2*dphi*(a+r*cos(theta)));
  return s;
}
-------------------------------------------------------------------------------------------------------------------------










PetscScalar dF4_3dBr(PetscScalar vpi1,PetscScalar vpi0,PetscScalar dr,PetscScalar vp,PetscScalar r)
{
  PetscScalar s;
  s = (vpi1-vpi0)/(2*dphi)+vp/r;
  return s;
}

PetscScalar dF4_3dBt(PetscScalar r,PetscScalar vpk1,PetscScalar vpk0,PetscScalar dphi)
{
  PetscScalar s;
  s = (vpk1-vpk0)/(2*r*dphi);
  return s;
}

PetscScalar dF4_3dBp(PetscScalar dt,PetscScalar vri1,PetscScalar vri0,PetscScalar dr,PetscScalar r,PetscScalar vr,PetscScalar vtk1,PetscScalar vtk0,PetscScalar dphi)
{
  PetscScalar s;
  s = -1/dt-(vri1-vri0)/(2*dr)-vr/r-(vtk1-vtk0)/(2*r*dphi);
  return s;
}

PetscScalar dF4_3dvr(PetscScalar Bpi1,PetscScalar Bpi0,PetscScalar dr,PetscScalar r,PetscScalar Bp)
{
  PetscScalar s;
  s = -(Bpi1-Bpi0)/(2*dr)-Bp/r;
  return s;
}

PetscScalar dF4_3dvt(PetscScalar r,PetscScalar Bpk1,PetscScalar Bpk0,PetscScalar dphi)
{
  PetscScalar s;
  s = -(Bpk1-Bpk0)/(2*dphi*r);
  return s;
}

PetscScalar dF4_3dvp(PetscScalar Bri1,PetscScalar Bri0,PetscScalar dr,PetscScalar Br,PetscScalar r,PetscScalar dphi,PetscScalar Btk1,PetscScalar Btk0)
{
  PetscScalar s;
  s = (Bri1-Bri0)/(2*dr)+Br/r+(Btk1-Btk0)/(2*dphi*r);
  return s;
}

PetscScalar dF4_3dBri1(PetscScalar vp,PetscScalar dr)
{
  PetscScalar s;
  s = vp/(2*dr);
  return s;
}

PetscScalar dF4_3dBri0(PetscScalar vp,PetscScalar dr)
{
  PetscScalar s;
  s = -vp/(2*dr);
  return s;
}

PetscScalar dF4_3dBtk1(PetscScalar r,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = vp/(r*2*dphi);
  return s;
}

PetscScalar dF4_3dBt0(PetscScalar r,PetscScalar vp,PetscScalar dphi)
{
  PetscScalar s;
  s = -vp/(r*2*dphi);
  return s;
}

PetscScalar dF4_3dBpi1(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = -vr/(2*dr);
  return s;
}

PetscScalar dF4_3dBpi0(PetscScalar vr,PetscScalar dr)
{
  PetscScalar s;
  s = vr/(2*dr);
  return s;
}

PetscScalar dF4_3dBpk1(PetscScalar r,PetscScalar vt,PetscScalar dphi)
{
  PetscScalar s;
  s = -vt/(2*r*dphi);
  return s;
}

PetscScalar dF4_3dBpk0(PetscScalar r,PetscScalar vt,PetscScalar dphi)
{
  PetscScalar s;
  s = vt/(2*r*dphi);
  return s;
}

PetscScalar dF4_3dvri1(PetscScalar Bp,PetscScalar dr)
{
  PetscScalar s;
  s = -Bp/(2*dr);
  return s;
}

PetscScalar dF4_3dvri0(PetscScalar Bp,PetscScalar dr)
{
  PetscScalar s;
  s = Bp/(2*dr);
  return s;
}

PetscScalar dF4_3dvtk1(PetscScalar Bp,PetscScalar dphi)
{
  PetscScalar s;
  s = -Bp/(2*dphi);
  return s;
}

PetscScalar dF4_3dvtk0(PetscScalar Bp,PetscScalar dphi)
{
  PetscScalar s;
  s = Bp/(2*dphi);
  return s;
}

PetscScalar dF4_3dvpi1(PetscScalar Br,PetscScalar dr)
{
  PetscScalar s;
  s = Br/(2*dr);
  return s;
}

PetscScalar dF4_3dvpi0(PetscScalar Br,PetscScalar dr)
{
  PetscScalar s;
  s = -Br/(2*dr);
  return s;
}

PetscScalar dF4_3dvpk1(PetscScalar r,PetscScalar Bt,PetscScalar dphi)
{
  PetscScalar s;
  s = Bt/(2*r*dphi);
  return s;
}

PetscScalar dF4_3dvpk0(PetscScalar r,PetscScalar Bt,PetscScalar dphi)
{
  PetscScalar s;
  s = -Bt/(2*r*dphi);
  return s;
}
---------------------------------------------------------------------------------------------------------------------------
