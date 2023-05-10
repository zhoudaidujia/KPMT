function [gama1,gama2]=gama(H00,H10,ita,E0)
[mx,my]=size(H00);
E=diag(ones(1,mx)*(E0+1i*ita));
I=diag(ones(1,mx));
t0=(E-H00)\(H10');
tt0=(E-H00)\H10;
gama1=t0;gama2=tt0;
tmpt1=t0;tmptt1=tt0;
tmpt3=I;tmptt3=I;
conver=1;conver2=1;
while conver>1.0d-7&&conver2>1.0d-7
    tmpt3=tmpt3*tmptt1;tmptt3=tmptt3*tmpt1;
    tmpII=(I-tmpt1*tmptt1-tmptt1*tmpt1);
    tmpt2=tmpII\(tmpt1*tmpt1);
    tmptt2=tmpII\(tmptt1*tmptt1);
    gama1=gama1+tmpt3*tmpt2;
    gama2=gama2+tmptt3*tmptt2;
    tmpt1=tmpt2;
    tmptt1=tmptt2;
    %----------------------------------------------------------------
	  conver=sum(sum(abs(tmpt2)));
	  conver2=sum(sum(abs(tmptt2)));
      %--------------------------------------------------------------
end
end
