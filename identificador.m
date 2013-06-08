source funciones.m;

[signal,fs,bps] = wavread('./grabaciones/all/8/69.wav');

%figure(3);
%plot(x);
%figure(4);
%plot(y);

inicio = 12e3;
n = 2*1024;
x = signal(inicio:inicio+n-1) .* hanning(n);

%como es un hombre tomar f0~200hz
muestras = floor(fs/200);
tx = liftrado(x,muestras);

%prediccion lineal
[a,err] = prediccion_lineal(x,14);
orden = 14;
a = a(:,orden);
freq = freqz( sqrt(err(orden)), [1;a],n,'whole');

df = fs/n;
f = [0:df:fs-df];


%por supuesto la culpa es de tu lsp
r = lsp2(a);
r2 =unique(abs(r))*fs/(2*pi);


figure(1);
clf;
%plot(f,abs(tx), 'k');
hold on;
plot(f,abs(freq), 'r');
stem(r2, 5*ones(size(r2)), 'k' );

figure(2);
plot(x);

r2'
