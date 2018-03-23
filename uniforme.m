% CUANTIFICACIÓN UNIFORME

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
%% 3) Cuantificación uniforme

%Parámetros cuantificador
N_niveles = 8;     %Número de niveles del cuantificador
M_din = 1;        %Márgen dinámico del cuantificador

partition = linspace(-M_din/2,M_din/2,N_niveles-1);
delta = partition(2)-partition(1);
codebook = -M_din/2-delta/2:delta:M_din/2+delta/2; 

[idx,xq] = quantiz(x1,partition,codebook); % Cuantifica
xq=xq';
figure 
plot(t1,x1,'x',t1,xq,'.')
legend('Original','Cuantificada','Location','NorthOutside');

%% 4) Prestaciones del cuantificador

% Valoración cualitativa de la cuantificación
if reproduce 
    player = audioplayer(xq,Fs);
    play(player);
end

% Valoración numérica a través del error cuadrático
e = x1 - xq;
figure
subplot(2,2,1), plot(t1,x1,t1,e);
title('Señal original vs error');
xlabel('tiempo (s)');
axis([6 12 -1 1]);

subplot(2,2,2), plot(t1,xq);
title('Señal cuantizada');
xlabel('tiempo (s)');
axis([6 12 -1 1]);

ec = e.^2;
subplot(2,2,3), plot(t1,ec);
title('Error cuadrático');
xlabel('tiempo (s)');
axis([6 12 0 10^(-2)]);

ecm = mean(ec);

% Ocurrencia de los símbolos
[prob, symbols] = hist(xq,unique(xq));
prob = prob./sum(prob);
prob = prob';
subplot(2,2,4), stem(symbols, 100.*prob);
title('Ocurrencia de los niveles de cuantificación');
xlabel('Nivel de cuantificación');
ylabel('Ocurrencia (%)');
axis([-M_din/2 M_din/2 0 50]);

% Tasa del cuantificador

%Tamaño a la entrada del cuantificador
tamx1 = length(x1)*16;              %16 bits por muestra

%El tamaño de la secuencia cuantificada sin comprimir es:
tamxq = length(xq)*ceil(log2(length(unique(xq))));

%La ganancia de compresión por cuantificación es:
gainq = tamx1/tamxq;              %Ganancia en la cuantificación

%% 5) Códificación y decodificación Huffman

%Aplicamos codificación Huffman a la señal cuantificada. 

% Construye el diccionario utilizando la ocurrencia de los símbolos:
[dict, avglen] = huffmandict(symbols,prob);
%dict es el diccionario
%avglen es la longitud media del código.
%Compruebe ambas variables.

% Codifica la secuencia utilizando el diccionario Hufmann que acabamos de
% fabricar
xcod = huffmanenco(xq, dict);

%Ahora podemos almacenar la secuencia, enviarla por un canal, etc...
ycod = xcod;  %Cuando estudiemos los canales, aquí sucederán calamidades...

% Decodifica la secuencia utilizando el diccionario
ydec = huffmandeco(ycod,dict);

% Comprueba que la secuencia decodificada en el receptor y la cuantificada
% son iguales
isequal(xq,ydec)

%% 6) Prestaciones en compresión Huffman

%El tamaño de la secuencia comprimida es
tamhuff1 = length(xcod);

%Aproximadamente lo mismo que:
tamhuff2 = length(xq)*avglen;

% Así que:
gainhuff1 = tamxq / tamhuff1;

%Podemos afirmar que la tasa de compresión conseguida es, en general:
gainhuff2 = log2(N_niveles)/avglen;

%% 8) Codificación y decodificación LZW
% Convertimos los niveles de cuantificación a símbolos binarios

mapping = [symbols [1:length(unique(xq))]'];
xqdec = xq;

for i = 1:length(unique(xq))
    I = xq == mapping(i,1);
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

%% 9) Prestaciones compresor LZW

%Tamaño de la secuencia recibida
tamlzw = length(ylzw)*13; %Las palabras código LZW tienen 13bits

%Ganancia de compresión
gainlzw = tamxq/tamlzw;

%% 10) Ganancias totales

gaintotal1 = gainq * gainhuff1;
gaintotal2 = gainq * gainlzw;
