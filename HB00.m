function [H]=HB00(r,a0,factor,B)
%Hk
        %area=((X2(J,1)-X2(I,1))*(X2(J,2)+X2(I,2))/2)*1.D-20;
        %exp(1i*factor*area*B);
%[rll,rl,rr,rrr,rc,rs]=coordinate(a0,a,N,M,0);
[mx,my]=size(r);
H=zeros(mx,mx);
Vppi=-2.66;
for j1=1:mx
    for j2=1:mx
        %----------------------------------------------
            tmpdx=r(j1,1)-r(j2,1);
            tmpdy=r(j1,2)-r(j2,2);
            tmpd=sqrt((tmpdx)^2+(tmpdy)^2);
            area=((r(j2,1)-r(j1,1))*(r(j2,2)+r(j1,2)));%*1.D-20;
            if tmpd<sqrt(3)*a0/3+0.1 && tmpd>0.01
                %H(j1,j2)=Vppi*exp(1i*factor*area*B);  
                H(j1,j2)=Vppi*exp(2i*pi/sqrt(3)/(a0^2)*B*area);
            elseif tmpd<0.01
                H(j1,j2)=0;
            end
    end
end 
end









