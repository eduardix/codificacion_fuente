% CUANTIFICACIÓN VECTORIAL

%% 1) Carga el archivo y realiza procesado previo
clear all
close all
clc

% Cambiar estos flags si no queremos reproducir o representar

representa = 1;  
reproduce = 0;

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
%% 3) Cuantificación y codificación vectorial
% Este cuantificador vectorial es de dos dimensiones. Es posible hacerlo de
% tantas dimensiones como se quiera. Con dos dimensiones es más sencillo
% apreciar cómo funciona

% Parámetros del cuantificador
dim = 2;        % Dimensiones (el presente código sólo opera con dim = 2)
offset = 1;     % Distancia entre dos muestras 
N_centr = 128;   % Número de centroides

[muestras,~]    = size(x1);
N_Parejas = muestras / dim;
X = zeros(N_Parejas,dim);
posicion = 1;

for count = 1:dim:muestras-offset    
    X(posicion,:) = [x1(count),x1(count+offset)];
    posicion = posicion + 1;    
end

%Cuantificación y codificación
[idx,ctrs] = kmeans(X,N_centr, 'Distance','sqEuclidean','start','cluster');
%idx = índices de los centroides. Cada centroide representa dos muestras.
%ctrs = valor de los centroides. 


%% 4) Representación del espacio de fases y los clusters
if representa 
    
    % Representamos las parejas de muestras
    figure;
    plot(X(:,1),X(:,2),'.')
    axis([-1 1 -1 1])
    xlabel('dimensión 1: muestra n');
    ylabel(['dimensión 2: muestra n + ' num2str(offset)]);    
    title('Espacio de fases');
    
    % Representa clusters y centroides
    figure;
    
    for i = 1:N_centr
        scatter(X(idx==i,1),X(idx==i,2),12)
        hold on   
    end

    if representa
        plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize',12,'LineWidth',2)
        plot(ctrs(:,1),ctrs(:,2),'ko','MarkerSize',12,'LineWidth',2)
    end
    axis([-1 1 -1 1]);
    xlabel('dimensión 1: muestra n');
    ylabel(['dimensión 2: muestra n + ' num2str(offset)]);    
    title('Espacio de fases particionado');  
    
end

% Decodificación a partir de idx
N_Parejas                           = length(idx);
Par_Muestras_Decod                  = zeros(N_Parejas,2);
posicion                            = 1:N_Parejas;
Par_Muestras_Decod(posicion,:)      = ctrs(idx(posicion),:);
Par_Muestras_Decod                  = Par_Muestras_Decod';
xq                                  = Par_Muestras_Decod(:);

%% 5) Error cuadrático medio o distorsión

% Valoración cualitativa de la cuantificación
if reproduce 
    player = audioplayer(xq,Fs);
    play(player);
end

figure;
% Valoración numérica a través del error cuadrático
e = x1 - xq;
subplot(2,2,1), plot(t1,x1,t1,e);
title('Señal original vs error');
xlabel('tiempo (s)');

subplot(2,2,2), plot(t1,xq);
title('Señal cuantificada');
xlabel('tiempo (s)');

ec = e.^2;
subplot(2,2,3), plot(t1,ec);
title('Error cuadrático');
xlabel('tiempo (s)');

ecm = mean(ec);

% Ocurrencia de los símbolos
[prob, symbols] = hist(idx,unique(idx));
prob = prob./sum(prob);
prob = prob';
subplot(2,2,4), stem(symbols, 100.*prob);
title('Ocurrencia de los símbolos');
xlabel('Símbolos');
ylabel('Ocurrencia (%)');

% Tasa del cuantificador

%Tamaño a la entrada del cuantificador
tamx1 = length(x1)*16;              %16 bits por muestra

%El tamaño de la secuencia cuantificada:
tamxq = length(idx)*log2(length(unique(idx)));

%La ganancia de compresión por cuantificación es:
gainq = tamx1/tamxq;              %Ganancia en la cuantificación 

%% 6) Código Huffman

% Construye el diccionario utilizando la ocurrencia de los símbolos:
[dict, avglen] = huffmandict(symbols,prob);
%dict es el diccionario
%avglen es la longitud media del código.
%Compruebe ambas variables.

% Codifica la secuencia utilizando el diccionario Hufmann que acabamos de
% fabricar
xhuffman = huffmanenco(idx, dict);

%Ahora podemos almacenar la secuencia, enviarla por un canal, etc...

ycod = xhuffman;
%% 7) Decodificación Huffman y vectorial

% El decodificador tiene dos diccionarios, el diccionario Huffman y el
% diccionario de los centroides. El primero sirve para decodificar el
% código Huffman. El segundo sirve para asignar valores concretos a las
% etiquetas del cuantificador vectorial.

% Decodifica la secuencia utilizando el diccionario
yq = huffmandeco(ycod,dict);

% Comprueba que ambas secuencias son iguales
isequal(idx,yq)

% Decodificación vectorial
N_Parejas                           = length(yq);
Par_Muestras_Decod                  = zeros(N_Parejas,2);
posicion                            = 1:N_Parejas;
Par_Muestras_Decod(posicion,:)      = ctrs(yq(posicion),:);
Par_Muestras_Decod                  = Par_Muestras_Decod';
ydec                                = Par_Muestras_Decod(:);

isequal(xq,ydec)

%% 7) Prestaciones de compresión Huffman

%Tomamos el tamaño de ydec, que es la secuencia recibida y decodificada
%tanto por el decodificador Hufmann como por el cuantificador vectorial no
%uniforme, y la comparamos con el tamaño de la secuencia xhuffman:

tamhuffman = length(ycod);
tamydec = length(idx)*log2(length(unique(idx)));

gainhuff = tamydec/tamhuffman;

%% 8) Compresor LZW
% Binarizamos y serializamos la secuencia:
xbin = dec2bin(idx);
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
tamlzw = length(xlzw)*13; %Las palabras código LZW tienen 13bits

%Ganancia de compresión
gainlzw = tamxq/tamlzw;

%% 10) Ganancias totales

gaintotal1 = gainq * gainhuff;
gaintotal2 = gainq * gainlzw;