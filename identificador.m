source funciones.m;
source ./mfcc/melfcc.m
inicio = [14,12,12,12,13,20,12,12,20,7]'*1e3;

vocales = [800 1600; 600 2500; 300 2500; 600 800; 300 800];
debe_dar = [vocales(4,:); vocales(1,:); vocales(2,:); vocales(3,:); vocales(4,:); vocales(5,:); vocales(1,:); vocales(2,:); vocales(3,:); vocales(4,:)];



for digito=[8]
	for L=[1:6]
	valores = [];
	tracto = [];
	lsp_coef = [];
		archivo = sprintf('./grabaciones/all/%d/%d.wav',digito,L);
		[signal,fs,bps] = wavread(archivo);

		n = 1024;
		[frames,t] = ventaneo(signal, n, 2, hanning(n));

		for K=[1:size(frames)(2)]
			x = frames(:,K);
			%prediccion lineal
			[a,err] = prediccion_lineal(x,14);
			orden = 14;
			a = a(:,orden);

			%por supuesto la culpa es de tu lsp
			r = lsp2(a);
			r2 =unique(abs(r))*fs/(2*pi);

			valores = [valores;r2'];

			tracto = [tracto,a];
			lsp_coef = [lsp_coef,r2];
		end

		clf;
		subplot(2,1,2);
		hold on;
		for K=[2:4]
			plot( t,lsp_coef(K,:), sprintf('-*%d',K-1) );
		end
		subplot(2,1,1);
		plot(signal);

		pause;
	end
end

%df = fs/n;
%f = [0:df:fs-df];
%f = f(1:n/2);
%for digito=[0,4,9]
%
%	freq = freqz(1, [1;tracto(:,digito+1)],n,'whole');
%	freq = freq(1:n/2);
%
%	digito
%	figure(1);
%	clf;
%	%plot(f,abs(tx), 'k');
%	hold on;
%	plot(log(f),log(abs(freq)), 'r');
%	stem(log(lsp_coef(:,digito+1)), 2*ones(size(lsp_coef(:,digito+1))), 'k' );
%	stem(log(lsp_coef(:,digito+1)), -10*ones(size(lsp_coef(:,digito+1))), 'k' );
%	stem( log(debe_dar(digito+1,:)), 10*ones(size(debe_dar(digito+1,:))), 'b');
%	stem( log(debe_dar(digito+1,:)), -10*ones(size(debe_dar(digito+1,:))), 'b');
%	pause;
%end
