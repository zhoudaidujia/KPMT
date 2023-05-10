function [rll,rl,rr,rrr,rc,rs]=coordinate(a0,a,N,M,az)
%=============================================================
%rl is the coordinate of left lead
%rr is the coordinate of right lead
%rc is the coordinate of center
%rs is the slice of center
%az means choose zigzag(1) or armchair(0) 
%az=1; %zigzag
r90=[cos(pi/2),-sin(pi/2);sin(pi/2),cos(pi/2)];
 Wz=sqrt(3)*a0/2;
rl=zeros(2*N,2);
rll=rl;
rr=zeros(2*N,2);
rrr=rr;
rc=zeros(2*N*M,2); %this is sample
rs=zeros(2*N,2,M); %this is slice of sample
rall=zeros(2*N*(M+4),2);
if az==1
    M=M+4;
y1=-1/sqrt(3)*a0-Wz*floor(N/2);y2=-1/sqrt(3)*a0+Wz*ceil(N/2);
x1=-0.25*a0-(M/2)*a0;x2=-0.25*a0+(M/2)*a0;
flag=0;
for j1=-500:1:500
    for j2=-500:1:500
        Ca=j1*a(1,:)+j2*a(2,:);
        Cb=Ca+(a(1,:)+2*a(2,:))/3;
        if Ca(1)<x2&&Ca(1)>x1&&Ca(2)>y1&&Ca(2)<y2
            flag=flag+1;
            rall(flag,:)=Ca;
        end
        if Cb(1)<x2&&Cb(1)>x1&&Cb(2)>y1&&Cb(2)<y2
            flag=flag+1;
            rall(flag,:)=Cb;
        end
    end
end
%=================================================================
rall=sortrows(roundn(rall,-5));
%=================================================================
elseif az==0
    M=M+4;
    Wa=a0/2;La=sqrt(3)*a0;
y1=-1/sqrt(3)*a0-M/2*La;y2=-1/sqrt(3)*a0+M/2*La;
x1=-0.25*a0-floor(N/2)*Wa;x2=-0.25*a0+ceil(N/2)*Wa;
flag=0;
for j1=-100:1:100
    for j2=-100:1:100
        Ca=j1*a(1,:)+j2*a(2,:);
        Cb=Ca+(a(1,:)+2*a(2,:))/3;
        if Ca(1)<x2&&Ca(1)>x1&&Ca(2)>y1&&Ca(2)<y2
            flag=flag+1;
            rall(flag,:)=Ca;
        end
        if Cb(1)<x2&&Cb(1)>x1&&Cb(2)>y1&&Cb(2)<y2
            flag=flag+1;
            rall(flag,:)=Cb;
        end
    end
end
rall=(r90*(rall'))';
%=======================================================
rall=sortrows(roundn(rall,-5));
%=======================================================
end
%rl=zeros(2*N,2);rll=rl;rr=zeros(2*N,2);rrr=rr;rc=zeros(2*N*M,2); %this is samplers=zeros(2*N,2,M); %this is slice of sample
rrr=rall(2*N*(M-1)+1:2*N*(M),:);
rall(2*N*(M-1)+1:2*N*(M),:)=[];
rr=rall(2*N*(M-2)+1:2*N*(M-1),:);
rall(2*N*(M-2)+1:2*N*(M-1),:)=[];
rll=rall(1:2*N,:);
rall(1:2*N,:)=[];
rl=rall(1:2*N,:);
rall(1:2*N,:)=[];
rc=rall;
for j=1:M-4;
    rs(1:2*N,:,j)=rall(1:2*N,:);
    rall(1:2*N,:)=[];
end

end
