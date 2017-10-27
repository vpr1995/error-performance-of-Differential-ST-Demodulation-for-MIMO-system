 
% Given data
M=8;
N=34;
fc = 800e6; % carrier frequncy in Hz
ns=1e5;     % number of samples
v = 22.352; % speed of vehicle in m/s
BW = 30e3;  % Band Width in Hz
c = 3e8;    % speed of light in m/s
n1 = c/fc;  % wavelength
maxfd = v/n1; %max doppler frequency in Hz
Ts = 1/BW;
 
sum1=0;sum2=0;sum3=0;sum4=0;
 
% Generating Channel Coefficients by using modified jakes simulator
h=hadamard(8); % hadamard matrix
for k=1:ns     % Number of channel samples
    m1=0;m2=0;m3=0;m4=0;
    for n=1:M  % Number of Oscillators
        on=(2*pi*n)/N;
        bn=(pi*n)/(M+1);
          m1 = m1 + (h(1,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*2)*bn));% for rayleigh channel 1 
          m2 = m2 + (h(2,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*3)*bn));% for rayleigh channel 2
          m3 = m3 + (h(3,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*4)*bn));% for rayleigh channel 3
          m4 = m4 + (h(4,n)*exp(1i*bn)*cos(2*pi*k*maxfd*Ts*cos(on)+(2*5)*bn));% for rayleigh channel 4
    end
    
    h1a(k)= m1;
    h2a(k)= m2;
    h3a(k)= m3;
    h4a(k)= m4;
 
    h1(k)= abs(m1);
    h2(k)= abs(m2);
    h3(k)= abs(m3);
    h4(k)= abs(m4);
  
    sum1=sum1+(h1(k)*h1(k));   %Power of channel 1
    sum2=sum2+(h2(k)*h2(k));   %Power of channel 2
    sum3=sum3+(h3(k)*h3(k));   %Power of channel 3
    sum4=sum4+(h4(k)*h4(k));   %Power of channel 4
 end
final1 = sum1/ns; %Average Power of channel 1
final2 = sum2/ns; %Average Power of channel 2
final3 = sum3/ns; %Average Power of channel 3
final4 = sum4/ns; %Average Power of channel 4

hp1= h1/sqrt(final1); %Normlized Channel Coefficients of channel 1
hp2= h2/sqrt(final2); %Normlized Channel Coefficients of channel 2
hp3= h3/sqrt(final3); %Normlized Channel Coefficients of channel 3
hp4= h4/sqrt(final4); %Normlized Channel Coefficients of channel 4
 
hp1c= h1a/sqrt(final1); %Normlized Channel Coefficients of channel 1
hp2c= h2a/sqrt(final2); %Normlized Channel Coefficients of channel 2
hp3c= h3a/sqrt(final3); %Normlized Channel Coefficients of channel 3
hp4c= h4a/sqrt(final4); %Normlized Channel Coefficients of channel 4
 
% plotting the channels 
t=1:400;
figure(1)
plot(t, 10*log(hp1(t)),'g');% rayleigh channel 1
hold on;
plot(t, 10*log(hp2(t)),'r');% rayleigh channel 2
hold on;
plot(t, 10*log(hp3(t)),'blue');% rayleigh channel 3
hold on;
plot(t, 10*log(hp4(t)),'black');% rayleigh channel 4
title('first 400 channel samples for the four rayleigh channels');
xlabel('Channel Samples');
ylabel('|h(t)| in dB'); 
legend('Channel 1','Channel 2','Channel 3','Channel 4');
 
% ***********************************************************************************************************************************************************
%  SER for Mr=1 case
snr = -5:10;
snrratio=10.^(snr./10);
Mt = 2;
Mr = 1;
iterations=1e5;

for i= 1:length(snr)
    s_old= sqrt(Mt)*eye(Mt);
    h1t = [hp1c(1);hp2c(1)];
    N1 = sqrt(1/2)*randn(Mt,Mr)+1i*sqrt(1/2)*randn(Mt,Mr);
    y1t_old = (sqrt(snrratio(i)/Mt).*s_old*h1t)+N1;
    err1=0;
    r1=randi(4);
    for l = 1:iterations
        codeword = sqrt(Mt)*[exp(1i*(r1-1)*pi/2) 0;0 exp(1i*(r1-1)*pi/2)];% codeword generation
        C=codeword;
        s_new = (1/sqrt(Mt))*(C*s_old);
        N1 = sqrt(1/2)*randn(Mt,Mr)+1i*sqrt(1/2)*randn(Mt,Mr);
        h1t = [hp1c(l); hp2c(l)];
        y1t_new = (sqrt(snrratio(i)/Mt)*(s_new*h1t))+N1;% transmitted signal
        n=zeros(1,4);
        for u=1:4% decoding
            codeword = sqrt(Mt)*[exp(1i*(u-1)*pi/2) 0;0 exp(1i*(u-1)*pi/2)];
            diff=y1t_new-((1/sqrt(Mt))*codeword*y1t_old);
            n(u)=(norm(diff,'fro'))^2;
        end
        [q,r]=min(n);
        if(r~=r1)
        err1=err1+1;
        end
    y1t_old=y1t_new;
    s_old=s_new;
    end
ser1(i)=err1/(iterations);% SER for Mr=1 case
end

figure(2)
semilogy(snr,ser1,'r');
title('Symbol error rate vs SNR for Mt=2 and Mr=1')
xlabel('SNR in db')
ylabel('Symbol Error Rate');
grid on;
 
% ***********************************************************************************************************************************************************
% % SER for Mr=2
snr = -5:10;
snrratio=10.^(snr./10);
Mt = 2;
Mr = 2;
iterations=1e5;


for i= 1:length(snr)
    s_old= sqrt(Mt)*eye(Mt);
    h1t = [hp1c(1) hp3c(1); hp2c(1) hp4c(1)];
    N2 = sqrt(1/2)*randn(Mt,Mr)+1i*sqrt(1/2)*randn(Mt,Mr);
    y1t_old = (sqrt(snrratio(i)/Mt).*s_old*h1t)+N2;
    err2=0;
    r1=randi(4);
    for l = 1:iterations
        codeword = sqrt(Mt)*[exp(1i*(r1-1)*pi/2) 0;0 exp(1i*(r1-1)*pi/2)];% codeword generation
        C=codeword;
        s_new = (1/sqrt(Mt))*(C*s_old);
        N2 = sqrt(1/2)*randn(Mt,Mr)+1i*sqrt(1/2)*randn(Mt,Mr);
        h1t = [hp1c(l) hp3c(l); hp2c(l) hp4c(l)];
        y1t_new = (sqrt(snrratio(i)/Mt)*(s_new*h1t))+N2;% transmitted signal
        n=zeros(1,4);
        for u=1:4% decoding
            codeword = sqrt(Mt)*[exp(1i*(u-1)*pi/2) 0;0 exp(1i*(u-1)*pi/2)];
            diff=y1t_new-((1/sqrt(Mt))*codeword*y1t_old);
            n(u)=(norm(diff,'fro'))^2;
        end
   [q,r]=min(n);
   if(r~=r1)
     err2=err2+1;
   end
   y1t_old=y1t_new;
   s_old=s_new;
end
ser2(i)=err2/(iterations); % SER for Mr=1 case
end 
 
figure(3)
semilogy(snr,ser2,'r');
title('Symbol error rate vs SNR for Mt=2 and Mr=1')
xlabel('SNR in db')
ylabel('Symbol Error Rate');
grid on;
 
figure(4)
semilogy(snr,ser1,'r',snr,ser2,'g');
title('Symbol error rate vs SNR for Mt=2 and Mr= 1 and 2')
legend('Mr=1','Mr=2');
xlabel('SNR in db')
ylabel('Symbol Error Rate');
grid on;