source funciones.m;
source ./mfcc/melfcc.m

rango = load('rango.txt');
patron = load('media.txt');
for digito=[0:9]
for L=[7]
	lsp_coef = [];
	archivo = sprintf('./grabaciones/all/%d/%d.wav',digito,L);
	[signal,fs,bps] = wavread(archivo);
	signal = signal(:,1);

	umbral = 0.5*dot(signal,signal)/length(signal);

	n = 1024;
	[frames,t] = ventaneo(signal, n, 2, hanning(n));

	for K=[1:size(frames)(2)]
		x = frames(:,K);
		energia = dot(x,x)/(n*0.374); %compensar el ventaneo

		%eliminar las ventanas de silencio
		if energia>umbral
			%prediccion lineal
			[a,err] = prediccion_lineal(x,14);
			orden = 14;
			a = a(:,orden);

			%por supuesto la culpa es de tu lsp
			r = lsp2(a);
			r2 =unique(abs(r))*fs/(2*pi);

			%si no tiene f0 no es vocal
			if( r2(2) < 800 )
				lsp_coef = [lsp_coef,r2(2:5)];
			end
		end
	end

	%tomar la parte central
	[rows,cols] = size(lsp_coef);
	lsp_coef = lsp_coef(:, floor(cols/4):cols-ceil(cols/4));
	%clasificar
	result = [];
	for K=[1:size(lsp_coef)(2)]
		%result = [result, clasificar(lsp_coef(:,K), rango)];
		result = [result, clasificar2(lsp_coef(:,K), patron)];
	end

	digito
	result
end
end

