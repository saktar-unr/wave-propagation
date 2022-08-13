N=500;                       %# samples 1D 
lambda= 0.5e-6;                 %wavelength
k= (2*pi)/lambda;                %wavenumber
f=1;
L=sqrt(N*lambda*f);                         %spatial grid side length
dx=L/N;                      %sample interval
x=linspace(-L/2, L/2-dx, N); %linear coordinates
[X, Y]= meshgrid(x,x);       %2D cordinates


z1=0;
A=ones(N);
U=A*exp(-1i*k*z1); %plane wave field 
theta=3.1416;
% U=A*exp(-1i*k*Y*sin(theta));
n=1.5;
d01=15e-3;
d02=15e-3;
h0=exp(-1i*k*n*(d01+d02));
t_lens=h0.*exp(-1i*k/(2*f)*(X.^2+Y.^2)); %lens
U2=U.*t_lens;

u2=propTF(U2,L,lambda,0.9);
I_b=abs(u2).^2;

lz=lambda*z;
I2=(w^2/lz)^2.*(jinc(w/lz*sqrt(X.^2+Y.^2))).^2;

figure
imagesc(x,x,I_b);
axis xy; axis square; colorbar; %display intensity
xlabel('x (m)'); ylabel('y (m)');

figure
plot(x,I_b(N/2+1,:)); xlabel('x (m)'); ylabel('I');

figure
imagesc(x,x,I2);
axis xy; axis square; colorbar; %display intensity
xlabel('x (m)'); ylabel('y (m)');

figure
plot(x,I2(N/2+1,:)); xlabel('x (m)'); ylabel('I');