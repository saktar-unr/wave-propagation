
L1=0.5;
M=250;
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;

lambda=0.5*10^-6; %wavelength
k=2*pi/lambda; %wavenumber
w=0.051; %source half widtsh (m)
z=2000; %propagation dist (m)

[X1,Y1]=meshgrid(x1,y1);
u1=rect(X1/(2*w)).*rect(Y1/(2*w)); %src field
I1=abs(u1.^2); %src irradiance

figure(1)
imagesc(x1,y1,I1);
axis square; set(gca,'YDir','normal')
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('z= 0 m');



function[out]=rect(x);
    %
    % rectangle function
    %
    % evaluates rect(x)
    % note: returns odd number of samples for full width
    %
    out=abs(x)<=1/2;
    end


u2=propTF(u1,L1,lambda,z); %propagation
x2=x1;
y2=y1;
I2=abs(u2.^2);

figure(2)
imagesc(x2,y2,I2);
axis square; set(gca,'YDir','normal')
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['z= ',num2str(z),' m']);

figure(3) %irradiance profile
plot(x2,I2(M/2+1,:));
xlabel('x (m)'); ylabel('Irradiance');
title(['z= ',num2str(z),' m']);

figure(4) %plot obs field mag
plot(x2,abs(u2(M/2+1,:)));
xlabel('x (m)'); ylabel('Magnitude');
title(['z= ',num2str(z),' m']);

figure(5) %plot obs field phase
plot(x2,unwrap(angle(u2(M/2+1,:))));
xlabel('x (m)'); ylabel('Phase (rad)');