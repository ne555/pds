1;

source dp.m
source dp2.m
source dpfast.m
source istft.m
source pvsample.m
source simmx.m

function r = autocorrelacion(y)
	r = autocor(y);
	%r = zeros(size(y));
	%n = length(r);
	%y = [y;zeros(n,1)];
	%for K=[1:n]
	%	r(K) = dot( y(1:n), y(K:n+K-1) );
	%end
end

function [R, r] = matrix_wiener_hopf(s,orden)
	r = autocorrelacion(s);
	%n = length(r)-1;
	n=orden;
	R = zeros(n,n);
	for K=[1:n]
		for L=[1:n]
			R(K,L) = r(abs(K-L)+1);
		end
	end

	r = -r(2:n+1);
end

function [a,G] = wiener_hopf(s,orden)
	[R,r] = matrix_wiener_hopf(s,orden);
	a = R\r;
	G = sqrt(r(1));
end

function [a, err] = prediccion_lineal(s,orden_maximo)
	%orden < length(s)
	r = autocorrelacion(s);
	a = zeros(orden_maximo, orden_maximo);
	err = zeros(orden_maximo,1);

	alpha = -r(2)/r(1);
	a(1,1) = alpha;
	err(1) = r(1) * (1-alpha**2);

	for K=[2:orden_maximo]
		alpha = -1/err(K-1) * ( r(K+1)+ dot(a([1:K-1],K-1),r([K:-1:2])) );
		a(K,K) = alpha;
		for L=[1:K-1]
			a(L,K) = a(L,K-1) + alpha*a(K-L,K-1);
		end
		err(K) = err(K-1)*(1-alpha**2);
	end

	%normalizar el error
	err/=r(1);
end

function [orden,graph] = error_prediccion_final(error_normalizado, umbral)
	graph = zeros(length(error_normalizado)-1,1);
	for K=[1:length(graph)]
		graph(K) = 1-error_normalizado(K+1)/error_normalizado(K);
	end

	orden = find(graph<umbral)(1);

end

function [orden,graph] = akaike(error_normalizado, n)
	Ip = zeros(size(error_normalizado));
	for K=[1:length(Ip)]
		Ip(K) = log(error_normalizado(K)) + 2*K/n;
	end
	graph=Ip;

	for orden=[2:length(Ip)-1]
		if( Ip(orden)<Ip(orden-1) && Ip(orden)<Ip(orden+1) )
			break;
		end
	end
end

function t = muestreo(inicial, final, frecuencia)
	dt = 1/frecuencia;
	t = [inicial:dt:final-dt]';
end


function [x,t] = senoidal(Amplitud, frecuencia, fase, ti, tf, fm)
	t = muestreo(ti,tf,fm);
	x = Amplitud*sin(2*pi*frecuencia*t + fase);
end

function [x,t] = sync(Amplitud, frecuencia, fase, ti, tf, fm)
	t = muestreo(ti,tf,fm);
	x = Amplitud*sinc(2*pi*frecuencia*t+fase);
end

function x = square_(value)
	x = 1*value>=0 - 1*(value<0);
end

function [x,t] = square(Amplitud, frecuencia, fase, ti, tf, fm)
	t = muestreo(ti,tf,fm);
	x = square_(aux);
end

function x = triangular_(t)
	if t>0.5
		x = 2*(1-t);
	else
		x = 2*t;
	end
	x -= .5;
	x *= 2;
end

function [x,t] = triangular(Amplitud, frecuencia, fase, ti, tf, fm)
	t = muestreo(ti,tf,fm);
	periodo = 1/frecuencia;
	x = zeros(size(t));
	for K = [1:length(t)]
		x(K) = Amplitud*triangular_( mod(t(K), periodo)/periodo  );
	end
	%x = Amplitud*mod(t+fase,periodo)-;
end

function [x,t] = dyrac(Amplitud, fase, ti, tf, fm)
	t = muestreo(ti,tf,fm);
	x = Amplitud*((t+fase)==0);
end

function x = ruido(Amplitud, ti, tf, fm)
	muestras = (tf-ti)*fm;
	x = Amplitud*(2*rand(muestras,1)-1);
end

function [h,theta] = respuesta_en_frecuencia(ma,ar,steps)
	theta = [0:pi/(steps-1):pi];
	z = exp(-j*theta);
	h = polyval(ma,z) ./ polyval(ar,z);
end

function c = cepstrum(x)
	fx = fft(x);
	c = ifft( log(abs(fx)) );
end

function fx = cep2four(c)
	fx = exp( fft(c) );
end

%devuelve en frecuencias
function [respuesta, excitacion] = liftrado(x, n)
	c = real( cepstrum(x) );
	s = length(c);

	%respuesta del tracto
	respuesta = zeros(s,1);
	respuesta(1:n) = c(1:n);
	respuesta(s-n+2:s) = c(s-n+2:s);

	%excitacion
	excitacion = zeros(s,1);
	excitacion(n+1:s-n+1) = c(n+1:s-n+1);

	respuesta = cep2four(respuesta);
	excitacion = cep2four(excitacion);
end

function count = zero_crossing(x)
	count = 0;
	for K=[1:length(x)-1]
		if( x(K)*x(K+1) < 0 )
			count++;
		end;
	end
end

function espectrograma(signal, tam, solap, fm)
%function espectrograma(signal, tam_vent, solap, fm)
%El solapamiento indica la cantidad de veces que se toma en cuenta una muestra (en las distintas transformadas)
	n = length(signal);
	fn = floor(tam/2)+1;
	step = floor(tam/solap);

	deltaf = fm/tam;

	%rellenar con 0 para la ultima ventana
	signal = [signal; zeros(mod(n,tam),1)];
	n = length(signal);
	s = [];
	for K=[1:step:n-tam]
		x = signal(K:K+tam-1); %extraccion
		%x .*= hamming(tam); %ventaneo
		x .*= gaussian(tam,.1);
		fx = fft(x)(1:fn); %como es simetrica solo me interesa una mitad
		s = [s abs(fx)];
	end

	f = [0;fm/2];
	t = [0;n-tam-1]/fm;

	imagesc( t, f, log(s) ); 
	%colormap('gray');
	%imagesc( [0:n], [0:deltaf:fm/2-deltaf], log(s) ); 
	set(gca,'ydir','normal');
	
end

function [P,Q] = lsp(ar)
	P = [1;ar+flipdim(ar);1];
	Q = [1;ar-flipdim(ar);-1];
end

function raiz = lsp2(ar)
	[P,Q] = lsp(ar);
	raiz = [roots(P);roots(Q)];
	%despreciar los conjugados
	[s,index] = sort(angle(raiz));
	index = index(15:30);

	%considerar solo los angulos
	%el modulo es 1
	raiz = angle(raiz(index));
end


function [result, time] = ventaneo(signal, tam, solap, ventana)
	n = length(signal);
	step = floor(tam/solap);

	%rellenar con 0 para la ultima ventana
	signal = [signal; zeros(mod(n,tam),1)];
	n = length(signal);
	result = [];
	for K=[1:step:n-tam]
		x = signal(K:K+tam-1); %extraccion
		x .*= ventana;
		result = [result x];
	end
	time = 1:step:n-tam;
end

function y = entre(a,b,c)
	if a<=b && b<=c
		y=1;
	else
		y=0;
	end
end

function voc = clasificar(x, rango)
	voc = -1;
	n = 5;
	for K=[1:n]
		if (entre(rango(K,1), x(1), rango(K,2)) || entre(rango(K+1,1), x(1), rango(K+1,2))) && (entre(rango(K,3), x(2), rango(K,3)) || entre(rango(K+1,3), x(2), rango(K+1,4))) && (entre(rango(K,5), x(3), rango(K,6)) || entre(rango(K+1,5), x(3), rango(K+1,6))) && (entre(rango(K,7), x(4), rango(K,8)) || entre(rango(K+1,7), x(4), rango(K+1,8)))
			%voc = floor((K-1)/2);
			voc = K;
		end
	end
end

function [frames2,t2] = cortar(signal)
	%[signal,fs,bps] = wavread(archivo);
	n = 1024;
	plot(signal);
	[frames,t] = ventaneo(signal, n, 2, hanning(n));
	a = max(abs(signal));
	m = (a^2)*0.2;
	frames2 = [];
	t2 = [];
	for i=1:size(frames)(2)
		aux=sum(frames(:,i).^2);
		if (aux>m)
			frames2 = [frames2 frames(:,i)]; 
			t2 = [t2 t(i)];
		end
	end
	
end

function signal2 = cortar2(signal)
	a = max(abs(signal));
	m = (a^2)*0.2;
	signal2 = [];
	for i=1:size(signal)
		if (abs(signal(i))>m)
			signal2 = [signal2; signal(i)]; 
		end
	end
	
end

function dif = diferencia(a, b)
	n = length(a);
	dif = 0;
	for i= 1:n
		dif = dif + abs(a(i)-b(i));
	end
	dif = dif/n;
end


function dif = comparar(m1,m2,m3,m4, num)
	%max = [2.05813857 12.0806228 15.0959306 245.95556];
	max = [1.697430	8.279271 10.857321 124.645135];
	m = [m1/max(1),m2/max(2),m3/max(3),m4/max(4)]';
	%m = [m1,m2,m3,m4]';
	distancia = [];
	for K=[0:4]
		%archivo = sprintf('./momentos/%d/%d.txt',num,K);
		archivo = sprintf('./patrones_nuevos512/patron%d-%d.txt',num,K);
		patron = load(archivo);
		for i = 1:4
			patron(:,i) = patron(:,i)/max(i);
		end
		SM = simmx(m,patron');
		[p,q,C] = dp(1-SM);
		distancia = [distancia;C(size(C,1),size(C,2))];
	end
	dif = min(distancia);
end

function [clase,costo] = clasificador(m1,m2,m3,m4, vocal)
	parecidos = inf*ones(10,1);

	contra = [1 6 6; 2 7 7; 3 8 8; 0 4 9; 5 5 5];
	for K=[1:3]
		parecidos( contra(vocal,K)+1 ) = comparar(m1,m2,m3,m4, contra(vocal,K) );
	end
	
	minimo = min(parecidos);
	clase=11;
	for i=1:10
		if (parecidos(i)==minimo)
			clase = i-1;
			break;
		end
	end
	costo = minimo;
end


function clasificado = clasificador_estadistico(signal, vocal, digito)
		umbral = 0.5*dot(signal,signal)/length(signal);		
		n = 512;
		[frames,t] = ventaneo(signal, n, 2, hanning(n));
		t2 = [];
		momentos = [];
		%transformada = [];
		for K=[1:size(frames)(2)]
			e = [];
			x = frames(:,K);
			energia = dot(x,x)/(n*0.374); %compensar el ventaneo
			%eliminar las ventanas de silencio
			if energia>umbral
				x = abs(fft(x));
				e = statistics(x);
				momentos = [momentos,e];
				%transformada = [transformada,x(1:length(x)/2)];
				t2 = [t2, t(K)];
			end
		end
		c = ceil(size(momentos)(2)/4);
		m = momentos(6:9,1:c)';
		%figure(2);
		%plot(transformada);
		%pause;
		
		%Versión real
		%[clasificado,costo] = clasificador(m(:,1),m(:,2),m(:,3),m(:,4), vocal); %Sin normalizar
	%	if vocal==5
	%		[clase2,costo2] = clasificador(m(:,1),m(:,2),m(:,3),m(:,4), 4); %Sin normalizar
	%		if costo2<costo
	%			clasificado = clase2;
	%		end
	%	end
		%clasificado = clasificador(m(:,1)/max(m(:,1)),m(:,2)/max(m(:,2)),m(:,3)/max(m(:,3)),m(:,4)/max(m(:,4)), vocal); %Normalizado
		
		
	%sprintf('Era: %d - Fue clasificado como: %d', digito, clasificado)
		% le decimos q vocal es así vemos q tal clasifica la consonante
		if (digito== 0 || digito == 4 || digito == 9)
				clasificado = clasificador(m(:,1),m(:,2),m(:,3),m(:,4), 4); %Sin normalizar
		end
		if (digito== 2 || digito == 7)
				clasificado = clasificador(m(:,1),m(:,2),m(:,3),m(:,4), 2); %Sin normalizar
		end
		if (digito== 1 || digito == 6)
				clasificado = clasificador(m(:,1),m(:,2),m(:,3),m(:,4), 1); %Sin normalizar
		end
		if (digito== 3 || digito == 8)
				clasificado = clasificador(m(:,1),m(:,2),m(:,3),m(:,4), 3); %Sin normalizar
		end
		if (digito == 5)
				clasificado = clasificador(m(:,1),m(:,2),m(:,3),m(:,4), 5); %Sin normalizar
		end

end

function voc = clasificar2(x, patron)
	distancia = [];
	for K=[1:size(patron)(1)]
		distancia = [distancia;norm(x-patron(K,:)')];
	end
	[descartar,voc] = min(distancia);
end


function y = mapear(x, mapa)
	y = zeros(size(x));
	for K=[1:length(x)]
		y(K) = mapa(x(K));
	end
end

function vocal = analizar( x )
	%mapa = [1 1 2 2 3 3 4 4 4 5 5];
	mapa =  [4 1 2 3 4 5 1 2 3 4];
	y = mapear(x,mapa);
	vocal = mode(y);
end

function vocal = dame_la_vocal( signal,fs,cant )
	patron = load('../base_datos/lsp.txt');
	[patron, factor] = normalizar(patron);
	lsp_coef = [];
	signal = signal(:,1); %no estereo

	umbral = 0.5*dot(signal,signal)/length(signal);
	n = 1024;
	[frames,t] = ventaneo(signal, n, 2, hanning(n));
	lsp_coef = caracteristicas_vocal(frames, fs, umbral, cant);

	%tomar la parte central
	[rows,cols] = size(lsp_coef);
	lsp_coef = lsp_coef(:, floor(cols/4):cols-ceil(cols/4));
	%clasificar
	result = [];
	for K=[1:size(lsp_coef)(2)]
		%result = [result, clasificar(lsp_coef(:,K), rango)];
		caracter = lsp_coef(:,K);
		caracter = caracter ./ factor';
		result = [result, clasificar2(caracter, patron)];
	end

	vocal = analizar(result);
end

function [mm, factor] = normalizar(m)
	factor = max(m);
	mm = zeros(size(m));
	for K=1:size(m)(2)
		mm(:,K) = m(:,K)/factor(K);
	end
end


function lsp_coef = caracteristicas_vocal(frames, fs, umbral, cant)
	lsp_coef = [];
	n = size(frames)(1);
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
			%r2 = r(2:cant+1)*fs/(2*pi);
			r2 = r(3:cant+2)*fs/(2*pi);
			f0 = r(1);

			%si no tiene f0 no es vocal
			%if( r2(2) < 800 )
			if( f0 < 800 )
				lsp_coef = [lsp_coef,r2];
			end
		end
	end
end

function [senal, fs] = extraer_senal(path, digito, n)
	pattern = cat(2, path, '%d/%d.wav'); 
	archivo = sprintf(pattern ,digito,n);
	[senal,fs,bps] = wavread(archivo);
	senal = senal(:,1);
end

function m = promedio_ponderado(senal)
	[rows,cols] = size(senal);
	senal = senal(:, floor(cols/4):cols-ceil(cols/4));

	m = zeros(1, size(senal)(1));
	n = size(senal)(2);

	%recta
	factor = linspace(0.5, 1.5, n);
	for K=[1:length(m)]
		for L=[1:n]
			m(K) += factor(L)*senal(K,L);
		end
	end
	m /= sum(factor);

end

function y = agregar_ruido(x, snr)
	n = length(x);
	ps = dot(x,x)/n;
	ruido = rand(size(x))-0.5;
	%normalizar el ruido
	pr = dot(ruido,ruido)/n;
	ruido /= sqrt(pr);
	%ajustar el ruido a snr
	pr = ps/ 10**(snr/10);

	ruido *= sqrt(pr);
	y = x+ruido;
end
