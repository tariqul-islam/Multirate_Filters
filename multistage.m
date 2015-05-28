clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%Design of Narrowband Filters using Multirate Technique
%
%
%This code takes an example
%
%First Filters it with a LPF filter
%
%Then Filter it with Three Stage Multirate Signal
%
%1006071 – Mohammad Tariqul Islam
%1006075 – Shuvo Newaz
%1006085 - Md. Maksudur Rahman
%1006087 – Md. Rafiul Amin
%1006090 - Ahnaf Sakib 
%1006190 - Moh Sabbir Saadat
%
%%%%%%%%%%%%%%%%%%%%%%%

%F1 = [1 10 50 100 300 1000 2000]; %frequencies in the signal
F1 = [1 10 50 100 300 1000 2000]; %frequencies in the signal
Fs = 8000; %Sampling Frequency
W = 2*pi*F1/Fs; %omega

%THE FUNCTION
f1 = @(n) (sin(W(1)*n+pi/4));
f2 = @(n) (sin(W(2)*n+2*pi/7));
f3 = @(n) (sin(W(3)*n+pi/6));
f4 = @(n) (sin(W(4)*n));
f5 = @(n) (sin(W(5)*n+pi/2));
f6 = @(n) (sin(W(6)*n+pi/3));
f7 = @(n) (sin(W(7)*n+7*pi/4));

ford = @(dp,ds,fp,fs,Fs) (ceil((-10*log10(dp*ds)-13)/(14.6*(fs-fp)/Fs)));

f = @(n) (f1(n)+f2(n)+f3(n)+f4(n)+f5(n)+f6(n)+f7(n));


delp = 10^-2;
dels = 10^-4;
delta = [delp dels];
fp = 75;
fs = 80;
fb = [fp fs];
m=[1 0];

fout = @(n) (f1(n)+f2(n)+f3(n));

n=0:1:2^16;
x = f(n);
true_filtered_output = fout(n);


%No Downsample/Upsample
[N, fo, mo, w] = firpmord(fb,m,delta,Fs); %filter order
N = ford(delta(1),delta(2),fp,fs,Fs);
%[N, fo, mo, w] = remezord(2*pi*fb/Fs,m,delta,2*pi);
%NK = kaiser_length(fb,m,delta,Fs);
b = firpm (N, fo, mo, w); %filter co eff
fprintf('Using FIR Filter: ')
tic
out = conv(x,b);
out = out(1:length(x));
toc
%Three Stage
M = 50;
L = 50;

M1 = 5;
M2 = 5;
M3 = 2;
L1 = 2;
L2 = 5;
L3 = 5;

%Filter Design
%first filter
FS1 = Fs;
fsd1 = Fs/M1 - Fs/(2*M);
[N1, fo, mo, w] = firpmord([fp fsd1],m,[delta(1)/6 delta(2)],FS1);
N1 = ford(delta(1)/6,delta(2),fp,fsd1,FS1);
b1D = firpm (N1, fo, mo, w);


%second filter
FS2 = Fs/M1;
fsd2 = Fs/(M1*M2) - Fs/(2*M);
[N2, fo, mo, w] = firpmord([fp fsd2],m,[delta(1)/6 delta(2)],FS2);
N2 = ford(delta(1)/6,delta(2),fp,fsd2,FS2);
b2D = firpm (N2, fo, mo, w);


%third filter
FS3 = Fs/(M1*M2);
fsd3 = Fs/(M1*M2*M3) - Fs/(2*M);
[N3, fo, mo, w] = firpmord([fp fsd3],m,[delta(1)/6 delta(2)],FS3);
N3 = ford(delta(1)/6,delta(2),fp,fsd3,FS3);
b3D = firpm (N3, fo, mo, w);

fprintf('Using Multistage Filters: ');
tic
%first downsample
K=length(b1D);
bx = fliplr(b1D);
for i=1:M1:length(x)
    x1((i-1)/M1+1) = sum(x(max(1,i-K+1):i).*bx(max(1,K-i+1):K));
end

%second downsample
K=length(b2D);
bx = fliplr(b2D);
for i=1:M2:length(x1)
    x2((i-1)/M2+1) = sum(x1(max(1,i-K+1):i).*bx(max(1,K-i+1):K));
end

%third Downsample
K=length(b3D);
bx = fliplr(b3D);
for i=1:M3:length(x2)
    x3((i-1)/M3+1) = sum(x2(max(1,i-K+1):i).*bx(max(1,K-i+1):K));
end

%first Upsample
x4 = upsample(x3,L1);
out4 = conv(x4,L1*b3D);
out4 = out4(1:length(x4));

%2nd Upsample
x5 = upsample(out4,L2);
out5 = conv(x5,L2*b2D);
out5 = out5(1:length(x5));

%3rd Upsample
x6 = upsample(out5,L3);
out6 = conv(x6,L3*b1D);
out6 = out6(1:min(length(x6),length(x)));
toc

%plots
%Time Domain
figure(1),subplot(4,1,1),plot(x);
title('Function');
xlim([1 10000]);
figure(1),subplot(4,1,2),plot(true_filtered_output);
title('True Output');
xlim([1 10000]);
figure(1),subplot(4,1,3),plot(out);
title('FIR Filter output');
xlim([1 10000]);
figure(1),subplot(4,1,4),plot(out6);
title('Multistage Multirate System Output');
xlim([1 10000]);

%frequency domain
figure(2),subplot(3,1,1),plot(0:Fs/length(true_filtered_output):Fs-Fs/length(true_filtered_output),abs(fft(true_filtered_output)));
xlim([0 500]);
title('True Output');
figure(2),subplot(3,1,2),plot(0:Fs/length(out):Fs-Fs/length(out),abs(fft(out)));
xlim([0 500]);
title('No Stage Output');
figure(2),subplot(3,1,3),plot(0:Fs/length(out6):Fs-Fs/length(out6),abs(fft(out6)));
xlim([0 500]);
title('Multi Stage Output');

%plotting the filters
figure(3),subplot(4,1,1);
plot(0:Fs/length(b):Fs-Fs/length(b),abs(fft(b)));
xlim([0 Fs/2]);
title('No Stage Filter');
xlabel('Frequency');
ylabel('Magnitude');
figure(3),subplot(4,1,2),plot(0:FS1/length(b1D):FS1-FS1/length(b1D),abs(fft(b1D)));
hold on,plot(FS2,1,'go');
xlim([0 FS1/2]);
title('First Filter');
xlabel('Frequency');
ylabel('Magnitude');
figure(3),subplot(4,1,3),plot(0:FS2/length(b2D):FS2-FS2/length(b2D),abs(fft(b2D)));
xlim([0 FS2/2]);
hold on,plot(FS3,1,'go');
title('Second Filter');
xlabel('Frequency');
ylabel('Magnitude');
figure(3),subplot(4,1,4),plot(0:FS3/length(b3D):FS3-FS3/length(b3D),abs(fft(b3D)));
xlim([0 FS3/2]);
title('Third Filter');
xlabel('Frequency');
ylabel('Magnitude');

%output of each stage
figure(4),subplot(3,2,2),plot(0:FS1/length(out6):FS1-FS1/length(out6),abs(fft(out6)));
title('Third Stage Upsample');
xlim([0 FS1/2]);
figure(4),subplot(3,2,1),plot(0:FS2/length(x1):FS2-FS2/length(x1),abs(fft(x1)));
xlim([0 FS2/2]);
title('First Stage Downsample')
figure(4),subplot(3,2,4),plot(0:FS2/length(out5):FS2-FS2/length(out5),abs(fft(out5)));
title('Second Stage Upsample');
xlim([0 FS2/2]);
figure(4),subplot(3,2,3),plot(0:FS3/length(x2):FS3-FS3/length(x2),abs(fft(x2)));
title('Second Stage Downsample');
xlim([0 FS3/2]);
figure(4),subplot(3,2,6),plot(0:FS3/length(out4):FS3-FS3/length(out4),abs(fft(out4)));
title('First Stage Upsample');
xlim([0 FS3/2]);
figure(4),figure(4),subplot(3,2,5),plot(0:FS3/M3/length(x3):FS3/M3-FS3/M3/length(x3),abs(fft(x3)));
title('Third Stage Downsample');
xlim([0 FS3/(2*M3)]);
