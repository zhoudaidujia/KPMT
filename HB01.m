function [H]=HB01(r1,r2,a0,factor,B)
%Hk
%[rll,rl,rr,rrr,rc,rs]=coordinate(a0,a,N,M,0);
[mx,my]=size(r1);
H=zeros(mx,mx);
Vppi=-2.66;
for j1=1:mx
    for j2=1:mx
        %----------------------------------------------
            tmpdx=r1(j1,1)-r2(j2,1);
            tmpdy=r1(j1,2)-r2(j2,2);
            tmpd=sqrt((tmpdx)^2+(tmpdy)^2);
            area=((r2(j2,1)-r1(j1,1))*(r2(j2,2)+r1(j1,2)));%*1.D-20;
            if tmpd<sqrt(3)*a0/3+0.1 && tmpd>0.01
                %H(j1,j2)=Vppi*exp(1i*factor*area*B);      
                H(j1,j2)=Vppi*exp(2i*pi/sqrt(3)/(a0^2)*B*area);
            elseif tmpd<0.01
                H(j1,j2)=0;
            end
    end
end 
end









