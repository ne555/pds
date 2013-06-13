source funciones.m;

patron = load('media3.txt');
for digito=[5]
for L=[1:7]
	lsp_coef = [];
	archivo = sprintf('./grabaciones/all/%d/%d.wav',digito,L);
	[signal,fs,bps] = wavread(archivo);
	vocal = dame_la_vocal(signal, fs);

	digito
	vocal
end
end

