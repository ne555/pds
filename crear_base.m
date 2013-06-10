source funciones.m;
inicio = [14,12,12,12,13,20,12,12,20,7]'*1e3;

vocales = [800 1600; 600 2500; 300 2500; 600 800; 300 800]; 
debe_dar = [vocales(4,:); vocales(1,:); vocales(2,:); vocales(3,:); vocales(4,:); vocales(5,:); vocales(1,:); vocales(2,:); vocales(3,:); vocales(4,:)];



%for digito=[0:9]
for digito=[7]

	m1 = [];
	m2 = [];
	m3 = [];
	m4 = [];
for L=[0:4]
valores = [];
tracto = [];
lsp_coef = [];
	archivo = sprintf('../grabaciones/all/%d/%d.wav',digito,10*L+digito+1);
	[signal,fs,bps] = wavread(archivo);
	signal=cortar2(signal);	
	n = ceil(length(signal)/16);
	[frames,t] = ventaneo(signal, n, 2, hanning(n));
	for K=[1:size(frames)(2)]
		x = frames(:,K);
		x = abs(fft(x));
		e = statistics(x);
		%prediccion lineal
		[a,err] = prediccion_lineal(x,14);
		orden = 14;
		a = a(:,orden);

		%por supuesto la culpa es de tu lsp
		r = lsp2(a);
		r2 =unique(abs(r))*fs/(2*pi);

		valores = [valores;r2'];

		tracto = [tracto,a];
		lsp_coef = [lsp_coef,e];
	end
	clf;
	subplot(2,1,2);
	hold on;
	for K=[6:8]
		lsp_coef(K,:)'
		plot( t,lsp_coef(K,:), sprintf('-*%d', K-5));
	end
	size(lsp_coef(6,:))(2)
    m1 = [m1 lsp_coef(6,1:30)'];
	m2 = [m2 lsp_coef(7,1:30)'];
	m3 = [m3 lsp_coef(8,1:30)'];
	m4 = [m4 lsp_coef(9,1:30)'];
	subplot(2,1,1);
	plot(signal);
	lsp_imp=lsp_coef(6:9,:)';
	save("m17.txt", "m1");
	save("m27.txt", "m2");
	save("m37.txt", "m3");
	save("m47.txt", "m4");

	%pause;
end
end
