function [GR]=sGFR(ita,E0,GR0,rs,rr,a0,factor,B)
[rx,ry,rz]=size(rs);
[H1]=HB00(rs(:,:,rz),a0,factor,B);
[H2]=H01(rs(:,:,rz),rr,a0);
[mx,my]=size(H1);
E=diag(ones(1,mx)*(E0+1i*ita));
%GRtmp1=inv(E-H1-H2*GR0*(H2'));
GRtmp1=(E-H1-H2/GR0*(H2'));
for j1=rz-1:-1:2
[H1]=HB00(rs(:,:,j1),a0,factor,B);
[H2]=HB01(rs(:,:,j1),rs(:,:,j1+1),a0,factor,B);
    %GRtmp2=inv(E-H1-H2*GRtmp1*(H2'));
    GRtmp2=(E-H1-H2/GRtmp1*(H2'));
    GRtmp1=GRtmp2;
end
%GR=inv(GRtmp1);
GR=(GRtmp1);
end
