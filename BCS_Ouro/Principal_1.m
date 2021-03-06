close all hidden
close all force
clear 
clc

% Definição dos pontos estacionários (indicação da posição (em
% amostras) dentro do exeperimento

P.Tinicio = 0;
P.Tfim = 1000;
P.Ts = 10; % Tempo de amostragem de referência para a reamostragem.
P.nsim = (P.Tfim-P.Tinicio)/P.Ts+1;
[P,dados] = Parametros_BCS(P);

entrada = [40,0.95];
[P,dados] = linear_moderno(entrada,P,dados)

%% simulação modelo nao-linear

P.Tinicio = 0;
P.Tfim = 10000;
P.Ts = 10; % Tempo de amostragem de referência para a reamostragem.
P.nsim = (P.Tfim-P.Tinicio)/P.Ts+1;
t = (P.Tinicio:P.Ts:P.Tfim)';

n = "malha aberta falha"
dados.falha = zeros(P.nsim,6);
Entradas.Frequencia = 40*ones(P.nsim,1);
Entradas.Abertura_valvula_choke = 0.95*ones(P.nsim,1);
Entradas.Pressao_reservatorio = P.p(1)*ones(length(Entradas.Frequencia),1);
Entradas.Pressao_manifold = P.p(6)*ones(length(Entradas.Frequencia),1);

Entradas.Frequencia(100:round(end/3)) = 39;
Entradas.Abertura_valvula_choke(round(end/2):round(2*end/3)) = 0.98;
[Med,Saidas,dados] = M_BCS(P,Entradas,dados,n)

N = Med(:,1);
Pc = Med(:,2);
Qc = Med(:,4);

%% Simulação de Modelo linear
U = [Entradas.Frequencia-dados.Estacionario.frequencia,Entradas.Abertura_valvula_choke-dados.Estacionario.abertura_valvula_choke];

[Y] = lsim(dados.ss.modelo_linear,U,t);

N_linear = Y(:,1)+dados.Estacionario.nivel;
Pc_linear = Y(:,2)+dados.Estacionario.pressao_choke;
Qc_linear = Y(:,3)+ dados.Estacionario.vazao_choke;

%% Figuras

figure
plot(t,N)
title('Nivel')
xlabel('segundos')
ylabel('metros')
hold on
plot (t,N_linear,'r:','LineWidth',2);
legend('Nonlinear','Linear')

figure
plot(t,Pc)
title('Pressao da choke')
xlabel('segundos')
ylabel('Bar')
hold on
plot (t,Pc_linear,'r:','LineWidth',2);
legend('Nonlinear','Linear')

figure
plot(t,Qc)
title('Vazão da choke')
xlabel('segundos')
ylabel('m^3/h')
hold on
plot (t,Qc_linear,'r:','LineWidth',2);
legend('Nonlinear','Linear')
%% Geração dos observadores UIO 
dados.K = 20;

O1 = gerador_UIO(dados.ss.modelo_linear,dados.ss.E(:,1),dados)
O2 = gerador_UIO(dados.ss.modelo_linear,dados.ss.E(:,2),dados)

%% Simulação de modelo não linear com disturbio
Entradas.Frequencia = 40*ones(P.nsim,1);
Entradas.Abertura_valvula_choke = 0.95*ones(P.nsim,1);
U = [Entradas.Frequencia-dados.Estacionario.frequencia,Entradas.Abertura_valvula_choke-dados.Estacionario.abertura_valvula_choke];

dados.falha(round(0.5*end/4):round(1.5*end/4),1) = -0.2;
dados.falha(round(2.5*end/4):round(3.5*end/4),2) = -0.2;

figure
plot(t,dados.falha(:,1:2),':','LineWidth',2)
xlabel('segundos','FontSize',20)
ylabel('Magntidude','FontSize',20)
legend('d_1','d_2')

[Med,Saidas,dados] = M_BCS(P,Entradas,dados,n);
S_nlinear = [Med(:,1) Med(:,2) Med(:,4)];
Saida = S_nlinear;
Y_desvio = [Med(:,1)-dados.Estacionario.nivel,Med(:,2)-dados.Estacionario.pressao_choke,Med(:,4)-dados.Estacionario.vazao_choke];

%% Resultado do UIO completo
[Y]=lsim(O1.sszy,[U Y_desvio],t);
Saida_01 = [Y(:,1)+dados.Estacionario.nivel,Y(:,2)+dados.Estacionario.pressao_choke,Y(:,3)+dados.Estacionario.vazao_choke];
Res_01 = Saida - Saida_01;

[Y]=lsim(O2.sszy,[U Y_desvio],t);
Saida_02 = [Y(:,1)+dados.Estacionario.nivel,Y(:,2)+dados.Estacionario.pressao_choke,Y(:,3)+dados.Estacionario.vazao_choke];
Res_02 = Saida - Saida_02;

figure
plot(t,Res_01)

figure
plot(t,Res_02)