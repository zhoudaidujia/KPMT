function [H]=H00(r,a0)
%Hk
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
            if tmpd<sqrt(3)*a0/3+0.1 && tmpd>0.01
                H(j1,j2)=Vppi;             
            elseif tmpd<0.01
                H(j1,j2)=0;
            end
    end
end 
end









