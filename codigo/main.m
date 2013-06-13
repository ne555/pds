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
	clase = clasificador_estadistico(signal, vocal, digito);

end
end
