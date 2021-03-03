clc;
theta_tl2=[];        
time=[];

fid = fopen('Coord_equil_E=0.05.txt','r'); %File_Name here

for m=1:1:20
    C=fscanf(fid, '%f', [3, 1728]); %1728 & 1692
end

for m=1:1:80 %Run for 80 time stamps
    
theta=[]; %For doing a mean of the left-right angles
C = fscanf(fid, '%f', [3, 1728]);
C=transpose(C);
x=C(1:end,1);
y=C(1:end,2);
z=C(1:end,3);
M=(z>15)&(z<50);
z_new=z(M);
x_new=x(M);
y_new=y(M);
N=(x_new>20)&(x_new<160);
z_new=z_new(N);
x_new=x_new(N);
y_new=y_new(N);
%scatter(x_new,z_new);   
z_n=[];
x_n=[];
th=0.5;
s=1;
for i=1:1:size(z_new,1)
        xi=(x_new(i));
        zi=(z_new(i));
        x1=(x_new)-xi*ones(size(x_new,1),1);
        z1=(z_new)-zi*ones(size(x_new,1),1);
        rho=sum(exp(-(x1.^2+z1.^2)/(2*s*s))./((2*3.1415*s)^(1.5)));
        if (rho>0.2)
            z_n=[z_n;zi];
            x_n=[x_n;xi];
        end 
end     
x_newer=x_n;
z_newer=z_n;
%scatter(x_newer,z_newer);
k1=boundary(x_newer,z_newer,0.5);
x_newer=x_newer(k1);
z_newer=z_newer(k1);
%scatter(x_newer,z_newer);
A=(z_newer>(min(z_newer))*1.15);
x_newer=x_newer(A);
z_newer=z_newer(A); 
%scatter(x_newer,z_newer);
%scatter(x_newer,z_newer);

% Fitting an ELLIPSE!
xc=x_newer-mean(x_newer);
yc=z_newer-mean(z_newer);
M=[xc.^2,yc.^2,xc.*yc,xc,yc];

% Fit ellipse through (xc,yc)
P0 = zeros(5,1);
fun = @(P) norm(M*P-1)^2;
nonlcon = @(P) nlcon(P);
P = fmincon(fun,P0,[],[],[],[],[],[],nonlcon);
xi = linspace(20,160);
yi = linspace(10,50);
[XI,YI]=meshgrid(xi-mean(x_newer),yi-mean(z_newer));
M=[XI(:).^2,YI(:).^2,XI(:).*YI(:),XI(:),YI(:)];
z=reshape(M*P-1,size(XI));

%%Uncomment to get plot of the projection
%  close all
%  h1=plot(x_newer,z_newer,'r.');
%  hold on
%  h2=contour(xi,yi,z,[0 0],'b');
%  axis equal
%  xlabel('X')
%  ylabel('Z')
%  legend('data','fit')
%  title('Plot of the X-Z Projection');

a=P(1);b=P(2);c=P(3);d=P(4);e=P(5);
ymin=23-mean(z_newer); %You set the z-height at which the contact angle will be found here it is 23

%Get ellipse parameters and use them to get the slope
b1=(d+c*ymin);
a1=a;
c1=(b*ymin^2+e*ymin-1);
xmin=(-b1-sqrt(b1^2-4*a1*c1))/(2*a1);
xmax=(-b1+sqrt(b1^2-4*a1*c1))/(2*a1);
dydx=-(d+2*a*xmin)/(2*b*ymin+c*xmin+e);
dydx2=-(d+2*a*xmax)/(2*b*ymin+c*xmax+e);
t1=(180/3.1415926)*atan(dydx);
t2=(180/3.1415926)*atan(dydx2);

if t1>0
theta=[theta;abs(t1)];
else 
theta=[theta;(180-abs(t1))];
end

if t2<0
theta=[theta;abs(t2)];
else 
theta=[theta;(180-abs(t2))];
end

m1=mean(theta);
  theta_tl2=[theta_tl2;m1];
end

plot(theta_tl2); %Plot the avg. theta's
legend('298');
hold off;
check = std(theta_tl2(61:71));
%csvwrite('Angle_e=1.csv',theta_tl);
