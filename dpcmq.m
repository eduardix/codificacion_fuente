% CUANTIFICACIÓN DPCM (predicción lineal)

%% 1) Carga el archivo y realiza procesado previo
clear all
close all
clc

% Cambiar estos flags si no queremos reproducir o representar

representa = 1;  
reproduce = 1;

[x,x1,Fs,Fsold,t,t1] = prev_process('hal9000.wav');

%Reproduce el sonido
if reproduce
    player = audioplayer(x1,Fs);
    play(player);
end

%% 2) Representa formas de onda y espectrograma 
if representa    
    figure; 
    subplot(2,2,1), plot(t,x(:,1));
    title('Fragmento original'), xlabel('tiempo(s)'), ylabel('amplitud');
    axis([0 18 -1 1]);
    
    subplot(2,2,3), plot(t1,x1(:,1));
    title('Fragmento tras filtro paso bajo'), xlabel('tiempo(s)'), ylabel('amplitud');
    axis([0 18 -1 1]);
    
    especparam.anchura_ventana = 64;
    especparam.noverlap = 0;
    especparam.nfft = 256;
    
    [S, F, T, P] = spectrogram(x1,especparam.anchura_ventana,especparam.noverlap,especparam.nfft,Fs);

    subplot(2,2,[2 4]), surf(T,F,10*log10(P),'edgecolor','none'); axis tight;
    view(0,90)
    xlabel('Tiempo (s)'); ylabel('Frecuencia (Hz)');
    title('Espectrograma (SPL,dB)');
    colorbar('north')    
end

%% 3) Cuantificación y DPCM
%Cuantificamos y codificamos mediante DPCM. En este apartado también
%veremos el resultado de la decodificación

%Parámetros cuantificador
N_niveles = 32;     %Número de niveles del cuantificador
M_din = 0.5;         %Márgen dinámico del cuantificador

%Predictor DPCM
acx = xcorr(x1,'coeff');
I = find(acx == 1);
a = acx(I+1);

predictor = [0 a]; % y[k]=a*x(n-1) -> Predictor de orden uno. 

partition   = linspace(-M_din/2,M_din/2,N_niveles-1);
delta       = partition(2)-partition(1);
codebook    = -M_din/2-delta/2:delta:M_din/2+delta/2;

% Codificamos usando DPCM
xenc_dpcm = dpcmenco(x1,codebook,partition,predictor);

% Decodificamos para comprobar prestaciones 
xdec = dpcmdeco(xenc_dpcm,codebook,predictor);
xdec = xdec';

%Muestra las señales
figure
subplot(2,2,1), plot(t1,x1);
xlabel('tiempo(s)');
ylabel('valor de la señal');
title('señal original');
subplot(2,2,3), plot(t1,xenc_dpcm);
xlabel('tiempo(s)');
ylabel('nivel de cuantificación')
title('señal DPCM');

subplot(2,2,[2 4]), plot(t1,x1,'x',t1,xdec,'.');
legend('Original','Decodificada','Location','NorthOutside');
xlabel('tiempo(s)');
ylabel('valor de la señal');
%% 4) Error cuadrático medio o distorsión

% Medimos las prestaciones de DPCM.

% Valoración cualitativa de la cuantificación y codificación
if reproduce 
    player = audioplayer(xdec,Fs);
    play(player);
end

% Valoración numérica a través del error cuadrático
e = x1 - xdec;
figure
subplot(2,2,1), plot(t1,x1,t1,e);
title('Señal original vs error');
xlabel('tiempo (s)');
%axis([0 18 -1 1]);

subplot(2,2,2), plot(t1,xdec);
title('Señal DPCM decodificada');
xlabel('tiempo (s)');
%axis([0 18 -1 1]);

ec = e.^2;
subplot(2,2,3), plot(t1,ec);
title('Error cuadrático');
xlabel('tiempo (s)');

ecm = mean(ec);

% Ocurrencia de los símbolos
[prob, symbols] = hist(xenc_dpcm,unique(xenc_dpcm));
prob = prob./sum(prob);
prob = prob';
subplot(2,2,4), stem(symbols, 100.*prob);
title('Ocurrencia de los símbolos');
xlabel('Símbolo');
ylabel('Ocurrencia (%)');

% Tasa del cuantificador

%Tamaño a la entrada del cuantificador
tamx1 = length(x1)*16;              %16 bits por muestra

%El tamaño de la secuencia cuantificada:
tamxq = length(xenc_dpcm)*ceil(log2(length(unique(xenc_dpcm))));

%La ganancia de compresión por cuantificación es:
gainq = tamx1/tamxq;              %Ganancia en la cuantificación 

%% 5) Codificación Huffman

%Aplicamos codificación Huffman a la señal cuantificada y codificada
%mediante DPCM. Esto nos permite pasar disminuir aun más la longitud media
%del código:

% Construye el diccionario utilizando la ocurrencia de los símbolos:
[dict, avglen] = huffmandict(symbols,prob);
%dict es el diccionario
%avglen es la longitud media del código.
%Compruebe ambas variables.

% Codifica la secuencia utilizando el diccionario Hufmann que acabamos de
% fabricar
xcod = huffmanenco(xenc_dpcm, dict);

%Ahora podemos almacenar la secuencia, enviarla por un canal, etc...
ycod = xcod;

% Decodifica la secuencia utilizando el diccionario
yenc_dpcm = huffmandeco(ycod,dict);

% Comprueba que ambas secuencias son iguales
isequal(xenc_dpcm,yenc_dpcm)

% Decodifica usando DPCM
ydec = dpcmdeco(yenc_dpcm,codebook,predictor);
ydec = ydec';

%Comprueba que las secuencias son iguales
isequal(ydec,xdec)

%% 5.5hide) Prestaciones en compresión Huffman

%El tamaño de la secuencia comprimida es
tamhuff1 = length(xcod);

%Lo mismo que:
tamhuff2 = length(xenc_dpcm)*avglen;

% Así que:
gainhuff1 = tamxq / tamhuff1;

%% 6) Codificación LZW

% Reconvertimos los símbolos

mapping = [unique(xenc_dpcm)' [1:length(unique(xenc_dpcm))]'];
xqdec = xenc_dpcm;

for i = 1:length(unique(xenc_dpcm))
    I = xenc_dpcm == mapping(i,1);
    xqdec(I) = mapping(i,2);
end

% Binarizamos y serializamos la secuencia:
xbin = dec2bin(xqdec);
xbin = xbin(:);

[xlzw tablalzw_tx] = norm2lzw(uint8(xbin));

%Ahora podemos almacenar la secuencia, enviarla por un canal, etc.

ylzw = xlzw; %Cuando estudiemos los canales, aquí sucederán calamidades...

%Decodificamos LZW
[ybin tablalzw_rx] = lzw2norm(uint16(ylzw));

%Comprobamos que ambas secuencias son iguales
isequal(xbin, char(ybin)')

%% 7hide) Prestaciones compresor LZW

%Tamaño de la secuencia recibida
tamlzw = length(xlzw)*13; %Las palabras código LZW tienen 13bits

%Ganancia de compresión
gainlzw = tamxq/tamlzw;

%% 8hide) Ganancias cuantificación + compresión

gaintotal1 = gainq * gainhuff1;
gaintotal2 = gainq * gainlzw;