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
%Then Filter it with Polyphase Techniques
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
Nx = length(x);
true_filtered_output = fout(n);

M=50;

%No Downsample/Upsample
[N, fo, mo, w] = firpmord(fb,m,delta,Fs); %filter order
N = ford(delta(1),delta(2),fp,fs,Fs);
%m = ceil(N/M);
%N = M*m-1;
b = firpm (N, fo, mo, w); %filter co eff
fprintf('With FIR Filter: ');
tic
out = conv(x,b);
out = out(1:length(x));
toc
b = [b zeros(1,ceil(N/M)*M-length(b))];
x = [x zeros(1,ceil(Nx/M)*M-length(x))];
N=N+1;
MM = ceil(N/M);
Mx = ceil(Nx/M);
filt_bank = zeros(M,MM);
x_bank = zeros(M,Mx);
yX = zeros(1,length(x_bank));

LB = length(x_bank);
for i=1:M
    filt_bank(M-i+1,:) = downsample(b(i:end),M);
    x_bank(i,:) = downsample(x(i:end),M);
end

fprintf('With Sequential Programming: ');
tic
for i=1:M
    y1 = conv(x_bank(i,:),filt_bank(i,:));
    y1 = y1(1:LB);
    yX = yX+y1; 
end
toc
y = zeros(1,length(x_bank));
fprintf('With Parallel Programming: ');
tic
parfor i=1:M
    y1 = conv(x_bank(i,:),filt_bank(i,:));
    y1 = y1(1:LB);
    y = y+y1;
end
toc
FsM = Fs/M;

y_up = upsample(y,M);
outx = conv(y_up,M*b);

figure(3),plot(0:Fs/length(outx):Fs-Fs/length(outx),abs(fft(outx)))
xlim([0 Fs/2]);
title('Filter Output')
%plots
%Time Domain
figure(1),subplot(4,1,1),plot(x);
title('Function');
xlim([1 10000]);
figure(1),subplot(4,1,2),plot(true_filtered_output);
title('True Output');
xlim([1 10000]);
figure(1),subplot(4,1,3),plot(out);
title('FIR Filter Output');
xlim([1 10000]);
figure(1),subplot(4,1,4),plot(outx);
title('Polyphase Output');
xlim([1 10000]);

figure(2),freqz(filt_bank(1,:));
title('Frequency Response of one of the filters.')