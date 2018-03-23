% Esta función abre el archivo de sonido y lo somete a un filtrado que está
% fuera del temario de la asignatura. Si desean conocer los detalles del 
% filtrado, pueden consultarme en mi despacho. Eduardo del Arco.

function [y,y1,Fs,Fsold,t,t1] = prev_process(string)    
    %Abre el archivo de sonido, idéntico al de la anterior práctica
    [y, Fs, ~, ~] = wavread(string);

    %Carga y aplica un filtro paso bajo para eliminar ruido HF
    load('fpb.mat')
    y1 = filter(fpb, y);

    %Toma un solo canal
    canal       = 1;                        
        
    diezma      = 4;    
    Fsold       = Fs;  
    Fs          = Fs/diezma;
    fragmento   = 3e5/diezma:5e5/diezma+1;    %Ventana temporal
    y1          = downsample(y1,diezma);        %Reducción de la frecuencia
    y1          = y1(fragmento,canal);
    t1 = fragmento/Fs;
    t = 1/Fsold:1/Fsold:length(y)/Fsold; 
end