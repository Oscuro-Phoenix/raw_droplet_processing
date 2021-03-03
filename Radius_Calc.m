clc;
%clear all;

fid = fopen('Coord_actuate_equil_E=0.175.txt','r'); %File_Name here
R=[]; %Array to store radius values

%Center of the simulation surface/box
cx=200; 
cy=200;



for m=1:1:20
    C=fscanf(fid, '%f', [3, 1728]); %1728 & 1692
end

for m=1:1:100
%Scan the file (assumes 1728 water molecules) change it if it is different
C = fscanf(fid, '%f', [3, 1728]);
C=transpose(C);

x=C(1:end,1);
y=C(1:end,2);
z=C(1:end,3);

B=(y>20)&(y<200);
x=x(B);
y=y(B);
z=z(B);
A=(x>20)&(x<200);
x=x(A);
y=y(A);
z=z(A);
F=(z>10.5)&(z<22);
x=x(F);
y=y(F);

%Plot the filtered projection of the droplet on the substrate
plot(x,y);

%Get the boundary points
k1=boundary(x,y,0.5);
x_new=x(k1);
y_new=y(k1);
xc=x_new-mean(x_new);
yc=y_new-mean(y_new);
M=[xc.^2,yc.^2,xc.*yc,xc,yc];

% Fit ellipse through (xc,yc)
P0 = zeros(5,1);
fun = @(P) norm(M*P-1)^2;
nonlcon = @(P) nlcon(P);
P = fmincon(fun,P0,[],[],[],[],[],[],nonlcon);
xi = linspace(20,200);
yi = linspace(20,200);
[XI,YI]=meshgrid(xi-mean(x_new),yi-mean(y_new));
M=[XI(:).^2,YI(:).^2,XI(:).*YI(:),XI(:),YI(:)];
z=reshape(M*P-1,size(XI));
 close all
 h1=plot(x_new,y_new,'r.');
 hold on
 h2=contour(xi,yi,z,[0 0],'b');
 axis equal
 xlabel('x')
 ylabel('y')
 legend('data','fit')
 a=P(1);b=P(2);c=P(4);d=P(5);
 R1=sqrt((1+(c*c*0.25/a)+(d*d*0.25/b))/a);
 x_cen1=(-c/(2*a))+R1;
 x_cen2=(-c/(2*a))-R1;
 if m==1
     cx=(-c/(2*a));
 end 
 
R=[R;2*sqrt((1+(c*c*0.25/a)+(d*d*0.25/b))/b)];
end

plot(R);
title('Plot of the contact diam. dynamics');
xlabel('Time (in ps)');
ylabel('Diameter (in A)');
legend('CD vs t');
hold off;
%csvwrite('Dia_e=1_T=298.csv',R); %Uncomment to store R values in a csv
%file for plotting