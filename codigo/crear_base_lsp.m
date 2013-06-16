
function crear_base_lsp(cant)
	source funciones.m;
	caracteristicas = [];
	for digito=[0:9]
		lsp_coef = [];
		for L=[0:4]
			%[signal, fs] = extraer_senal('../grabaciones/entrenamiento/', digito, L);
			[signal, fs] = extraer_senal('../otro/', digito, 10*L+digito+1);

			umbral = 0.5*dot(signal,signal)/length(signal);

			n = 1024;
			[frames,t] = ventaneo(signal, n, 2, hanning(n));

			lsp_coef = [lsp_coef; promedio_ponderado(caracteristicas_vocal(frames, fs, umbral, cant))];
		end

		caracteristicas = [caracteristicas; mean(lsp_coef)];
	end

	save('-ascii', '../base_datos/lsp.txt', 'caracteristicas');
end
