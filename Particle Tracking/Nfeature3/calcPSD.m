% calculate PSD by Welch's method
Fs = 1000;   % sampling frequency
t = 0:1/Fs:.296;
x = 10*cos(2*pi*t*20)+randn(size(t));
plot(t,x)
pwelch(x,[],[],[],Fs);
[Pxx,w]=pwelch(x,[],[],[],Fs);
