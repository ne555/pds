source funciones.m;

patron = load('media.txt');
for digito=[0:9]
for L=[7]
	lsp_coef = [];
	archivo = sprintf('./grabaciones/all/%d/%d.wav',digito,L);
	[signal,fs,bps] = wavread(archivo);
	vocal = dame_la_vocal(signal, fs);

	digito
	vocal
end
end

