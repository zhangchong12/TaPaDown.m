%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%this is the main function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%=====================================================%%%%%%%%%%%
%%%to prepare the code
clc
clear all
%%%to load the data and rearrange the data
load gravity.csv
gra_air=gravity(:,3);
bou=reshape(gra_air,518,179);
bou1=flipud(bou);
bou1=rot90(bou1,-1);
dg2=bou1(:,1:179);
%%%if the data are not that smooth and with lots of noises
% dg2=wiener2(dg2,[2,2]);
%%%the grid numbers and grid intervals
dx=0.4;
dy=0.4;
[lx ly]=size(dg2);
x=0:1:lx-1;
y=0:1:ly-1;
[xx,yy]=meshgrid(0.4*x,0.4*y);
%%%the wavelegth and wavenumber
longx=dx*(lx-1);
longy=dy*(ly-1);
k_frequencytotal=frequencytotal(lx,ly,longx,longy);
%%%to prevent the zero value
k_frequencytotal(1,1)=1*10^(-6); 
%%%the observe level of the data 
z=0;
%%%the grid point interval of the downward continuation 
depth_n=5;
%%%the depth of the downward continuation 
z_up=-0.4*depth_n;  
% % % % % % %%%%%这里把输入的原始数据作为向下延拓的真值
% % % dg2_down=dg2;
% % % dg2=real(ifft2(fft2(dg2).*exp(-(z-z_up).*k_frequencytotal)));
%%%%%%%%%%%=====================================================%%%%%%%%%%%
%%%to calculate the vertical derivatives by the FFT 
for n1=1:4;
F_dg_z2=fft2(dg2).*(k_frequencytotal).^n1;
fft_dg_z2=real(ifft2(F_dg_z2));  
fft_dg_z2n(:,:,n1)=fft_dg_z2;
end
%%%to calculate the potential of the field by the FFT 
F_Ug2=fft2(dg2)./k_frequencytotal;
Ug2=real(ifft2(F_Ug2));%%%%%%%%%%
Ug2_temp=Ug2;
%%%to calculate the vertical derivatives by the ISVD
for n2=1:2:3
Ug2_kuox(1,:)=Ug2_temp(1,:);%%%%%%%%expand the x-direction
Ug2_kuox(2:lx+1,:)=Ug2_temp;
Ug2_kuox(lx+2,:)=Ug2_temp(lx,:);
Ug2_kuoy(:,1)=Ug2_temp(:,1);%%%%%%%%expand the y-direction
Ug2_kuoy(:,2:ly+1)=Ug2_temp;
Ug2_kuoy(:,ly+2)=Ug2_temp(:,ly);

chafen_U_xx2=(Ug2_kuox(3:lx+2,:)-2.*Ug2_kuox(2:lx+1,:)+Ug2_kuox(1:lx,:))./(dx).^2;
chafen_U_yy2=(Ug2_kuoy(:,3:ly+2)-2.*Ug2_kuoy(:,2:ly+1)+Ug2_kuoy(:,1:ly))./(dy).^2;

chafen_U_zz2=-(chafen_U_yy2+chafen_U_xx2);
Ug2_temp=chafen_U_zz2;
chafen_dg_z2n(:,:,n2)=chafen_U_zz2;
end

dg2_temp=dg2;
for n3=2:2:2
dg2_kuox(1,:)=dg2_temp(1,:);
dg2_kuox(2:lx+1,:)=dg2_temp;
dg2_kuox(lx+2,:)=dg2_temp(lx,:);
dg2_kuoy(:,1)=dg2_temp(:,1);
dg2_kuoy(:,2:ly+1)=dg2_temp;
dg2_kuoy(:,ly+2)=dg2_temp(:,ly);

chafen_dg2_xx2=(dg2_kuox(3:lx+2,:)-2.*dg2_kuox(2:lx+1,:)+dg2_kuox(1:lx,:))./(dx).^2;
chafen_dg2_yy2=(dg2_kuoy(:,3:ly+2)-2.*dg2_kuoy(:,2:ly+1)+dg2_kuoy(:,1:ly))./(dy).^2;

chafen_dg2_zz2=-(chafen_dg2_yy2+chafen_dg2_xx2);
dg2_temp=chafen_dg2_zz2;
chafen_dg2_zz2n(:,:,n3)=chafen_dg2_zz2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%#####################DOWNWARD CONTINUATION#####################%%%%%%%%%%%
%%%1、downward continuaiton by the FFT FFT method
TF_dg0=dg2+(z-z_up).^1.*fft_dg_z2n(:,:,1)./prod(1:1)+(z-z_up).^2.*fft_dg_z2n(:,:,2)./prod(1:2);

%%%1、downward continuaiton by the ISVD Taylor method
TISVD_dg0=dg2+(z-z_up).^1.*(chafen_dg_z2n(:,:,1)+chafen_dg2_zz2n(:,:,1))./prod(1:1)+(z-z_up).^2.*(chafen_dg_z2n(:,:,2)+chafen_dg2_zz2n(:,:,2))./prod(1:2);

%%%1、downward continuaiton by the ISVD Pade method
b0=dg2;
b1=chafen_dg_z2n(:,:,1)+chafen_dg2_zz2n(:,:,1);
b2=(chafen_dg_z2n(:,:,2)+chafen_dg2_zz2n(:,:,2))./prod(1:2);
bb=b0;
bb(bb==0)=eps;   
pp0=b0; 
qq0=1;
qq1=-b1./bb;
qq2=-(-b1.^2 + b0.*b2)./bb.^2;
PadeISVD_q=(qq0+(z-z_up).*qq1+(z-z_up).^2.*qq2);
PadeISVD_q(PadeISVD_q==0)=eps;
PadeISVD_dg0=(pp0)./PadeISVD_q; 
% %%%%%%%%%%%=====================================================%%%%%%%%%%%
save TF_dg0; 
save TISVD_dg0; 
save PadeISVD_dg0;
%%%%%%%%%%%=====================================================%%%%%%%%%%%
%%%MAPPING
figure(1)
subplot(2,2,1)
contourf(xx,yy,dg2);
k1=caxis;
title('(a) gravity anomalies of downward')
xlabel('x(m)')
ylabel('y(m)')
hco1=colorbar ;
t1=get(hco1,'YTickLabel');
t1=strcat(t1,'mGal');
set(hco1,'YTickLabel',t1);

subplot(2,2,2)
contourf(xx,yy,TF_dg0);
% caxis(k1);
title('(d) downward continuation by FFT-Taylor')
xlabel('x(m)')
ylabel('y(m)')
hco1=colorbar ;
t1=get(hco1,'YTickLabel');
t1=strcat(t1,'mGal');
set(hco1,'YTickLabel',t1);

subplot(2,2,3)
contourf(xx,yy,TISVD_dg0);
% caxis(k1);
title('(c) downward continuation by ISVD-Taylor')
xlabel('x(m)')
ylabel('y(m)')
hco1=colorbar ;
t1=get(hco1,'YTickLabel');
t1=strcat(t1,'mGal');
set(hco1,'YTickLabel',t1);

subplot(2,2,4)
contourf(xx,yy,PadeISVD_dg0);
% caxis(k1);
title('(d) downward continuation by ISVD-Pade')
xlabel('x(m)')
ylabel('y(m)')
hco1=colorbar ;
t1=get(hco1,'YTickLabel');
t1=strcat(t1,'mGal');
set(hco1,'YTickLabel',t1);