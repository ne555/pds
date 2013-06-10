source funciones.m;

for digito=[0:9]
%for digito=[7]
	m1 = [];
	m2 = [];
	m3 = [];
	m4 = [];
for L=[6:7]
	momentos = [];
	%archivo = sprintf('../grabaciones/all/%d/%d.wav',digito,10*L+digito+1);
	archivo = sprintf('./grabaciones/all/%d/%d.wav',digito,L);
	[signal,fs,bps] = wavread(archivo);
	vocal = dame_la_vocal(signal,fs);

	signal=cortar2(signal);	
	n = ceil(length(signal)/16);
	[frames,t] = ventaneo(signal, n, 2, hanning(n));
	for K=[1:size(frames)(2)]
		x = frames(:,K);
		x = abs(fft(x));
		e = statistics(x);
		momentos = [momentos,e];
	end
	clf;

    m1 = momentos(6,1:30)';
	m2 = momentos(7,1:30)';
	m3 = momentos(8,1:30)';
	m4 = momentos(9,1:30)';


	clasificado = clasificador(m1,m2,m3,m4, vocal);
	
	sprintf('Era: %d - Fue clasificado como: %d', digito, clasificado)
	%subplot(2,1,1);
	%plot(signal);
	
	%pause;
end
end
