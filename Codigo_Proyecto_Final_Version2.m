% Proyecto final Mecanica Vibratoria Grupo 7
clc
clear
%close all

Habilitar_Grafico_de_modos_de_vibracion = 0;
Habilitar_Grafico_de_respuestas_modales = 0;
Habilitar_Grafico_de_respuestas = 0;
Habilitar_Grafico_de_Movimiento_Relativo_tiempo = 1;
Habilitar_Grafico_de_Movimiento_Relativo_y_x = 0;
Habilitar_Grafico_de_Movimiento_Relativo_y_x_Secuencial = 0;
Habilitar_Movimiento_Theta = 0;
Calculo_y_grafico_trasmisibilidad = 0;


%--------Datos del problema-------
E = 70e9; #Modulo de young aluminio; 200e9; #Modulo de young SAE1045
MIy = 48.163e-9;  #Momento de inercia X
MIx = 12.305e-9;  #Momento de inercia Y
LH = 0.26;  #Longitud del perfil horizontal
LV = 0.35;  #Longitud del perfil vertical

#---Masas---
mh = 0.10634; #masa del perfil horizontal aluminio
mv = 0.28315; #masa del perfil vertical aluminio
masa_base = 4.15; #masa de la base
masa_Hotend = 0.1109; #masa del hotend

#---Rigideces---
k_perfil_v_X = (E*MIx)/LV;
k_perfil_v_Y = (E*MIy)/LV;
k_perfil_v_Torsional = 5.3e7;
##k_taco_X = 299.2;
##k_taco_Y = 299.8;
k_taco_X = 59.985;
k_taco_Y = 59.999;

#---Amortiguamiento---
c_perfil_v_X = 0;
c_perfil_v_Y = 0;
c_perfil_v_Torsional = 0;
c_taco_X = 10;
c_taco_Y = 10; %Estos C no son utilizados en los calculos, se hicieron los analisis de amortiguamiento variando dseta

#---Condiciones Iniciales---
  #Posiciones
Xo1 = 0.01;
Yo1 = 0.01;
Xo2 = 0;
Yo2 = 0;
Theta = 0;

  #Velocidades
Xo1p = 0.01;
Yo1p = 0;
Xo2p = 0;
Yo2p = 0;
Thetap = 0;

#---Tipos de vibraciones y cargas---
Vibraciones = 'Libres';    %Colocar 'Libres' o 'Forzadas'
Carga = 'Senoidal';        %Colocar 'Senoidal', 'Cosenoidal' o 'Exponencial'
AmplitudX = 0.005;  #Amplitud del movimiento de la base en X
AmplitudY = 0.005;
wx = 10;
wy = 10;


#---Ordenamiento de valores en vectores---

Xo = [Xo1;Xo2;Yo1;Yo2;Theta];          %Condiciones iniciales de posicion
Xop = [Xo1p;Xo2p;Yo1p;Yo2p;Thetap];            %Condiciones iniciales de velocidad

N = 1000;

%------Matrices y vectores auxiliares----
y = linspace(0 , 1, 5+1); #Sirve para graficar los phi normalizados
t = linspace(0 , 10*10 , N); #Si quiero que los graficos muestren menos tiempo, cambio el valor de Tp
taux = linspace(0 , 1*10 , N);  #Vector de tiempo para graficos de movimientos relativos Y(X)
I = eye(5);
vector_ceros = zeros(1,5);  %Esto lo uso para agregarlo a la matriz phi despues para considerar el empotramiento

%-------Matrices-------
#------#
m  = [masa_base,mv,mh+masa_Hotend];

M = [   m(1)   ,    0    ,    0    ,    0    ,    0    ;
         0     ,m(2)+m(3),    0    ,    0    ,    0    ;
         0     ,    0    ,   m(1)  ,    0    ,    0    ;
         0     ,    0    ,    0    ,m(2)+m(3), m(3)*LH ;
         0     ,    0    ,    0    , m(3)*LH ,m(3)*(LH^2)];

#------#
k  = [k_taco_X,k_taco_Y,k_perfil_v_X,k_perfil_v_Y,k_perfil_v_Torsional];

K = [k(1)+k(3) ,  -k(3)  ,      0      ,    0    ,    0;
      -k(3)    ,   k(3)  ,      0      ,    0    ,    0;
        0      ,    0    ,  k(2)+k(4)  ,  -k(4)  ,    0;
        0      ,    0    ,    -k(4)    ,   k(4)  ,    0;
        0      ,    0    ,      0      ,    0    ,   k(5)];

#------#

c  = [c_taco_X,c_taco_Y,c_perfil_v_X,c_perfil_v_Y,c_perfil_v_Torsional];

C = [c(1)+c(3) ,  -c(3)  ,      0      ,    0    ,    0;
      -c(3)    ,   c(3)  ,      0      ,    0    ,    0;
        0      ,    0    ,  c(2)+c(4)  ,  -c(4)  ,    0;
        0      ,    0    ,    -c(4)    ,   c(4)  ,    0;
        0      ,    0    ,      0      ,    0    ,   c(5)];  %Estas matriz fue previamente calculada en hoja con metodo de Lagrange



%---------Autovectores y Autovalores--------
[AutoV,lamda] = (eig(K,M));         %Calculo Autovectores normalizados con respecto a la masa y autovalores
Wn = diag(lamda.^0.5);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            %Calculo frecuencias naturales de cada elemento
%------Saco el phi normalizado a si mismo para ver los modos---------
for i = 1:length(AutoV(:,1))
  phi(:,i) = AutoV(:,i)/max(abs(AutoV(:,i)));
endfor

%--Calculo mediante metodo modal--

Wng = [((k_taco_X/m(1))^0.5),((k_perfil_v_X/m(2))^0.5),((k_taco_Y/m(1))^0.5),((k_perfil_v_Y/m(2))^0.5),((k_perfil_v_Torsional/m(3))^0.5)];
AutoVT = transpose(AutoV);    %Transpuesta del autovector
Mm = AutoVT*M*AutoV;  #Matriz masa modal
Km = AutoVT*K*AutoV;  #Matriz rigidez modal
Cm = AutoVT*C*AutoV;  #Matriz amortiguacion modal
Xom = Mm/I*AutoVT*M*Xo; #Vector CI de posicion modal
Xopm = Mm/I*AutoVT*M*Xop; #Vector CI velocodad modal

for i = 1:5
  dseta(i) = Cm(i,i)/(2*Wn(i)*Mm(i,i)); %Vector dseta calculado con matrices modales
endfor

Ax = AmplitudX*((k_taco_X^2)+(c_taco_X*wx)^2)^0.5;  #Fuerza transmitida de la base a la masa 1 en X
Ay = AmplitudY*((k_taco_Y^2)+(c_taco_Y*wy)^2)^0.5;  #Fuerza transmitida de la base a la masa 1 en Y

alfax = atan(-(c_taco_X*wx)/k_taco_X);
alfay = atan(-(c_taco_Y*wy)/k_taco_Y);

w = [wx;0;wy;0;0];  #vector frecuencias de movimiento de la base
A = [Ax;0;Ay;0;0];  #vector frecuencias de amplitudes de la base
alfa = [alfax;0;alfay;0;0]; #vector frecuencias de desfasajes de la base

%----------Calculos segun el tipo de carga definido previamente---------
switch Vibraciones
  case 'Forzadas'
    An = AutoVT*A;  %Modal de amplitudes
    alfan = AutoVT*alfa;  %Modal de alfas
    wn = AutoVT*w;  %Modal de frecuencias
    switch Carga
      case 'Senoidal'
        for i = 1:length(An)
          for j = 1: length(t)
            Xn(i,j) = (An(i)/((Km(i,i)-Mm(i,i)*(wn(i))^2)^2 + (Cm(i,i)*wn(i))^2)^0.5)*sin(wn(i)*t(j) - atan((Cm(i,i)*w(i))/(Km(i,i)-Mm(i,i)*(wn(i)^2))) - alfan(i));
          endfor
        endfor
        X = AutoV*Xn;

      case 'Cosenoidal'
        for i = 1:length(Pon)
          for j = 1: length(t)
            Xn(i,j) = (An(i)/((Km(i,i)-Mm(i,i)*(wn(i))^2)^2 + (Cm(i,i)*wn(i))^2)^0.5)*cos(wn(i)*t(j) - atan((Cm(i,i)*w(i))/(Km(i,i)-Mm(i,i)*(wn(i)^2))) - alfan(i));
          endfor
        endfor
        X = AutoV*Xn;
      otherwise
        fprintf('Error en el ingreso de Tipo de carga');
    endswitch
  case 'Libres'
    for i = 1:length(dseta)
      for j = 1:length(t)
        Xn(i,j) = (exp(-dseta(i)*Wn(i)*t(j)))*(Xom(i)*cos(((1-dseta(i)^2)^0.5)*Wn(i)*t(j))+((Xopm(i)+dseta(i)*Wn(i)*Xom(i))/(Wn(i)*(1-dseta(i)^2)^0.5))*sin(((1-dseta(i)^2)^0.5)*Wn(i)*t(j)));
      endfor
    endfor
    X = AutoV*Xn;

 otherwise
   fprintf('Error en el ingreso de Tipo de vibracion')
endswitch

#----Defino vectores para facilitar la interpretacion del codigo----
  #Posiciones absolutas
  X1 = X(1,:);
  X2 = X(2,:);
  Y1 = X(3,:);
  Y2 = X(4,:);
  Y3 = X(5,:)*LH + Y2; #Paso coordenadas rotacionales a traslacionales

  #Posiciones relativas entre el extrusor y la cama caliente
  X31 = X2-X1;
  Y31 = Y3-Y1;
  maxX31 = max(X31)
  maxY31 = max(Y31)


%--------Valores de trasmisibilidad modal-----------

if Calculo_y_grafico_trasmisibilidad == 1

TrasmisibilidadX = max(X1)/AmplitudX
TrasmisibilidadY = max(Y1)/AmplitudY

##  j=1;
##  for i = 0:0.01:length(t)
##    Tdx(j) = ((i^2 + (c_taco_X*Wng(1))^2)/((i-m(1)*Wng(1)^2)^2 + (c_taco_X*Wng(1))^2))^0.5;
##    rx(j) = abs(Wng(1)/((i/m(1))^0.5));
##    if j <= 10
##      Tdx(j)
##    endif
##    if j >= 3
##      if (Tdx(j-1) > Tdx(j-2) && Tdx(j-1) > Tdx(j))
##        Kmcriticox = i
##      endif
##    endif
##    j++;
##  endfor
##
##  j=1;
##  for i = 0:0.01:length(t)
##    Tdy(j) = ((i^2 + (Cm(3,3)*wn(3))^2)/((i-Mm(3,3)*wn(3)^2)^2 + (Cm(3,3)*wn(3))^2))^0.5;
##    ry(j) = abs(wn(3)/((i/Mm(3,3))^0.5));
##    if j >= 3
##      if (Tdy(j-1) > Tdy(j-2) && Tdy(j-1) > Tdy(j))
##        Kmcriticoy = i
##      endif
##    endif
##    j++;
##  endfor
##
  Kmcriticox = (wn(1)^2)*Mm(1,1);
  Kmcriticoy = (wn(3)^2)*Mm(3,3);
##
##  #----Calculo "K critico" para de los tacos----"
##
##
##      Ktrasmisibilidad = Km;
##      Ktrasmisibilidad(1,1) = Kmcriticox;
##      Ktrasmisibilidad(3,3) = Kmcriticoy;
##      Ktrasmisibilidadreal = inv(AutoVT)*Ktrasmisibilidad*inv(AutoV)
##      Kcriticox = Ktrasmisibilidadreal(1,1)-k(3)
##      Kcriticoy = Ktrasmisibilidadreal(3,3)-k(4)
##
##      Ktrasmisibilidad = Km;
##      Ktrasmisibilidad(3,3) = Kmcriticoy;
##      Ktrasmisibilidadreal = inv(AutoVT)*Ktrasmisibilidad*inv(AutoV);
##      Kcriticoy = Ktrasmisibilidadreal(3,3)-k(4);
##
##
##  fprintf('Rigidez critica del taco en la direccion de X = %.3f\n', Kcriticox);
##  fprintf('Rigidez critica del taco en la direccion de Y = %.3f\n', Kcriticoy);
##
##  figure(8);
##    subplot(2,1,1)
##    plot(rx,Tdx,'k');
##    xlabel('r = w/wn');
##    ylabel('Td = X/Y');
##    title('Trasmisibilidad');
##
##    subplot(2,1,2)
##    plot(ry,Tdy,'k');
##    xlabel('r = w/wn');
##    ylabel('Td = X/Y');
##    title('Trasmisibilidad');

endif
%--------Grafico los modos normalizados-------

phi = [vector_ceros; phi];    %Agrego fila de ceros a la matriz de phi para que tenga en cuenta la parte fija (empotramiento)

if Habilitar_Grafico_de_modos_de_vibracion == 1
  figure(1);

  subplot(1,5,1);
  plot(phi(:,1),y,'k');
  xlabel('Phi');
  ylabel('Elemento');
  title('Modo 1');

  subplot(1,5,2);
  plot(phi(:,2),y,'k');
  xlabel('Phi');
  ylabel('Elemento');
  title('Modo 2');

  subplot(1,5,3);
  plot(phi(:,3),y,'k');
  xlabel('Phi');
  ylabel('Elemento');
  title('Modo 3');

  subplot(1,5,4);
  plot(phi(:,4),y,'k');
  xlabel('Phi');
  ylabel('Elemento');
  title('Modo 4');

  subplot(1,5,5);
  plot(phi(:,5),y,'k');
  xlabel('Phi');
  ylabel('Elemento');
  title('Modo 5');

endif

%--------Grafico respuestas modales para confirmar que dieron ondas senoidales puras------

if Habilitar_Grafico_de_respuestas_modales == 1

  figure(2)

  subplot(5,1,1);
  plot(t,Xn(1,:), 'k');
  xlabel('tiempo');
  ylabel('1');
  title('Respuesta modal 1');

  subplot(5,1,2);
  plot(t,Xn(2,:), 'k');
  xlabel('tiempo');
  ylabel('2');
  title('Respuesta modal 2');


  subplot(5,1,3);
  plot(t,Xn(3,:), 'k');
  xlabel('tiempo');
  ylabel('3');
  title('Respuesta modal 3');

  subplot(5,1,4);
  plot(t,Xn(4,:), 'k');
  xlabel('tiempo');
  ylabel('4');
  title('Respuesta modal 4');

  subplot(5,1,5);
  plot(t,Xn(5,:), 'k');
  xlabel('tiempo');
  ylabel('5');
  title('Respuesta modal 5');

endif

%--------Grafico respuestas de las diferentes masas-------

if Habilitar_Grafico_de_respuestas == 1

  figure(3);

  subplot(5,1,1);
  plot(t,X1, 'b');
  xlabel('tiempo');
  ylabel('X1');
  title('Respuesta abs de la m1 en X');

  subplot(5,1,2);
  plot(t,Y1, 'b');
  xlabel('tiempo');
  ylabel('Y1');
  title('Respuesta abs de la m1 en Y');

  subplot(5,1,3);
  plot(t,X2, 'r');
  xlabel('tiempo');
  ylabel('X2');
  title('Respuesta abs de la m2 en X');

  subplot(5,1,4);
  plot(t,Y2, 'r');
  xlabel('tiempo');
  ylabel('Y2');
  title('Respuesta abs de la m2 en Y');

  subplot(5,1,5);
  plot(t,Y3, 'g');
  xlabel('tiempo');
  ylabel('Y3');
  title('Respuesta abs de la m3 en Y');

endif

#---- Grafico de movimientos relativos entre el extrusor y la cama caliente ----

if Habilitar_Grafico_de_Movimiento_Relativo_tiempo == 1

  figure(4);
  subplot(2,1,1);
  plot(taux,X31, 'b');
  xlabel('tiempo');
  ylabel('Mov. relativo X');
  title('Movimiento rel. entre m1 y m3 en X');

  subplot(2,1,2);
  plot(taux,Y31, 'r');
  xlabel('tiempo');
  ylabel('Mov. relativo Y');
  title('Movimiento rel. entre m1 y m3 en Y');

endif

#---- Grafico el movimiento relativo Y en funcion de X ----

for i = 1:(N/10)
  X31aux(i) = X31(i);
  Y31aux(i) = Y31(i);
endfor

if Habilitar_Grafico_de_Movimiento_Relativo_y_x == 1

  figure(5)
  plot(X31aux,Y31aux, 'k');
  xlabel('Mov. relativo X');
  ylabel('Mov. relativo Y');
  title('Movimiento relativo resultante');

endif


%---- Grafico secuencialmente el movimiento relativo Y en funcion de X, para esto tengo que habilitarlo previamente ----

if Habilitar_Grafico_de_Movimiento_Relativo_y_x_Secuencial == 1
  for i = 1:length(X31aux)

    figure(6)
    plot(X31aux(i),Y31aux(i), 'k');
    xlabel('Mov. relativo X');
    ylabel('Mov. relativo Y');
    title('Movimiento relativo resultante');

    pause(0.02); % Agregar una pausa para visualizar cada punto
    hold on;
  endfor
endif

if Habilitar_Movimiento_Theta == 1

  figure(7)
  plot(taux,X(5,:)*LH, 'k');
  xlabel('Mov. relativo X');
  ylabel('Mov. relativo Y');
  title('Movimiento relativo resultante');

endif
