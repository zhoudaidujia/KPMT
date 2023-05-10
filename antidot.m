clc;clear;
%-------------------------------------------------------------
%factor=-1.519*10^19;	%						!(2*Pi*e/h)=(2*3.1415926*1.602D-19)/6.626D-34
%-------------------------------------------------------------
e=1.602176487*10^-19;
h=6.62606896*10^-34;
factor=e/h*pi;
a0=2.49;
N=400; %width of ribbon even for zigzag
M=40; %the number of unit cells of ribbon even for all
ita=10^-4;
a=a0.*[1,0;
     -0.5,sqrt(3)/2];
 Natom=2*N*M;
Natomall=2*N*(M+4);
 Wz=sqrt(3)*N*a0/2; Wa=N*a0/2; Lz=M*a0; La=sqrt(3)*M*a0;
 tic
[rll,rl,rr,rrr,rc,rs]=coordinate(a0,a,N,M,1);
'nime1'
toc
%--------------------------------------------------------------------
[H1]=H00(rll,a0);
[H2]=H01(rll,rl,a0);
[mx,my]=size(H1);
ne=100;
nb=100;
Energy=linspace(-0.4,0,ne);
Con=zeros(1,ne);
Con2=zeros(nb,ne);
Bm=linspace(0,0.008,nb);
flagstep=0;
for jj=1:ne
    E0=Energy(jj);
    tic
    %[gama1,gama2]=gama(H1,H2,ita,E0);
    [gama1,gama2]=gama(H1,H2,ita,E0);
    %B=Bm(jjb);
    flagstepmag=0;
for jjb=1:nb
    time1=clock;
    %E0=Energy(jj);
    B=Bm(jjb);
%--------------------------------------------------------------------
%left and right lead surface Green's function
E=diag(ones(1,mx)*(E0+1i*ita));
%I=diag(ones(1,mx));
%GL=inv(E-H1-(H2')*gama2);
GL=(E-H1-(H2')*gama2);
%GR0=inv(E-H1-H2*gama1);
GR0=(E-H1-H2*gama1);
[GR]=sGFR(ita,E0,GR0,rs,rr,a0,factor,B);
%---------------------------------------------------------------------
%total Green's function
[Hh1]=HB00(rs(:,:,1),a0,factor,B);
[Hh01]=H01(rl,rs(:,:,1),a0);
[Hh12]=HB01(rs(:,:,1),rs(:,:,2),a0,factor,B);
SL=(Hh01')/GL*Hh01;
SR=Hh12/GR*(Hh12');
%G11=inv(E-Hh1-SL-SR);
G11=(E-Hh1-SL-SR);
GamaL=1i*(SL-SL');
GamaR=1i*(SR-SR');
%Con(jj)=real(2*e^2/h*trace(GamaL*G11*GamaR*(G11')));
Con2(jjb,jj)=real(trace(GamaL/G11*GamaR/(G11')));
Con(jj)=real(trace(GamaL/G11*GamaR/(G11')));
time2=clock;
timemag=etime(time2,time1)
flagstepmag=flagstepmag+1
end
    flagstep=flagstep+1
    timefinish=toc
end
s1='antidot_ENE';
name=s1;
save (name,'Con2')












