source funciones.m;

result = zeros(10,10);
lavocal = zeros(10,10);
for digito=[0:9]
	for L=[0:9]
		lsp_coef = [];
		%archivo = sprintf('../grabaciones/all/%d/%d.wav',digito,L);
		archivo = sprintf('../grabaciones/seleccionados/%d/%d.wav',digito,L);
		[signal,fs,bps] = wavread(archivo);
		signal = signal(:,1);
		signal = agregar_ruido(signal, 0);
		vocal = dame_la_vocal(signal,fs, 8);

		clase = clasificador_estadistico(signal,vocal, digito);

		lavocal(digito+1, L+1) = vocal;
		result(digito+1, L+1) = clase;
	end
end

confusion = zeros(10,10);
for K=[1:size(result)(1)]
	for L=[1:size(result)(2)]
		confusion(K,result(K,L)+1)++;
	end
end

