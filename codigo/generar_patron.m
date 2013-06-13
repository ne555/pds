source funciones.m;
for digito=[0:9]
for L=[1:6]
	lsp_coef = [];
	archivo = sprintf('./grabaciones/all/%d/%d.wav',digito,L);
	[signal,fs,bps] = wavread(archivo);


	signal=cortar2(signal);	
	%n = ceil(length(signal)/51);
	n = 512;
	momentos = [];
	[frames,t] = ventaneo(signal, n, 2, hanning(n));
	for K=[1:size(frames)(2)]
		x = frames(:,K);
		x = abs(fft(x));
		e = statistics(x);
		momentos = [momentos,e];
	end
	%m1 = momentos(6,1:100)';
	%m2 = momentos(7,1:100)';
	%m3 = momentos(8,1:100)';
	%m4 = momentos(9,1:100)';

	m = momentos(6:9,:)';
	salida = sprintf('./momentos/%d/%d.txt',digito,L);
	save('-ascii', salida, 'm');
end
end


