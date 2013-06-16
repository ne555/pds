source funciones.m;

result = zeros(10,15);

for digito=[0:9]
%for digito=[7]
	m1 = [];
	m2 = [];
	m3 = [];
	m4 = [];
	for L=[5:14]
		momentos = [];
		%archivo = sprintf('../grabaciones/all/%d/%d.wav',digito,10*L+digito+1);
		%archivo = sprintf('../grabaciones/prueba/%d/%d.wav',digito,L);
		[signal, fs] = extraer_senal('../otro/', digito, 10*L+digito+1);
		%signal = agregar_ruido(signal, 30);
		%[signal,fs,bps] = wavread(archivo);
		vocal = dame_la_vocal(signal,fs, 6);
		%clase = clasificador_estadistico(signal, vocal, digito);

		result(digito+1, L+1) = vocal;
	end
end
