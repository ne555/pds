1;

function r = autocorrelacion(y)
	r = zeros(size(y));
	n = length(r);
	y = [y;zeros(n,1)];
	for K=[1:n]
		r(K) = dot( y(1:n), y(K:n+K-1) );
	end
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
	raiz = [angle(roots(P));angle(roots(Q))];
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

function voc = clasificar2(x, patron)
	distancia = [];
	for K=[1:size(patron)(1)]
		distancia = [distancia;norm(x-patron(K,:)')];
	end
	[descartar,voc] = min(distancia);
end
