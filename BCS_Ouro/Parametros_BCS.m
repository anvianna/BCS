function [P,dados] = Parametros_BCS(P)

         %% Parâmetros de simulação
        dados = struct;
        %% Parâmetros de construção do sistema LEA
    
        P.rho = [855 855];                  % [kg/m^3]            Massa específica do fluido nos volumes de controle  [CV1 CV2]         
        P.h = [32 23.2 0];                  % [m]                 Pronfundidade do reservatório, bomba e choke    
        P.pm =  101325;                     % [Pa]                Pressão do manifold equivalente a uma atmosfera
        P.k_choke = 0.228065321489704;      % resultado da otimização
        P.k_resv  = 0.720849094248127;      % resultado da otimização
        P.g = 9.81;                         % [m/s^2]             Aceleração da gravidade
        P.PI =  0.85e-8;                    % [m^3/(Pa*s)]        Indice de produtividade do poço
        P.mu = [0.00888667 0.00888667];     % [Pa.s]              Viscosidade do fluido nos volumes de controle [CV1 CV2]
        P.B = [1.8000e+04 1.8000e+04];      % [Bar]                Compressibilidade do fluido [CV1 CV2]
        P.B0 = [0 0];                       % [adm] Constante de fricção B0
        P.B1 = [0 0];                       % [adm] Constante de fricção B1                
        P.r = [0.11 0.0375];                % [m]                 Raio dos volumes de controle [CV1 CV2]
        P.l = [9.3 22.7];                   % [m]                 Comprimentos de tubo entre reservatório e a bomba e entre a bomba e a choke
        P.V = [P.l(1)*pi*P.r(1)^2, ...
               P.l(2)*pi*P.r(2)^2];         % [m^3]               Volumes dos volumes de controle [CV1 CV2]
        P.Aanular = 0.033595;               % [m^2]               Área do Anular

        %% Valores médios

        P.l_bar = (P.l(1) + P.l(2))/2;                                    % [m]       Comprimento médio
        P.r_bar = ((P.r(1)*P.l(1))+(P.r(2)*P.l(2)))/(P.l(1)+P.l(2));      % [m]       Raio médio
        P.A_bar =  pi*P.r_bar^2;                                          % [m^2]     Área média
        P.rho_bar = ((P.rho(1)*P.V(1))+(P.rho(2)*P.V(2)))/(P.V(1)+P.V(2));% [kg/m^3]  Massa específica média
        P.med_Val = [P.l_bar P.r_bar P.A_bar P.rho_bar];                  %           Valores médios 

        %% Parâmetros da bomba e da válvula
        P.Head = [-86.085161541490240,10.078950606511444,1.843640642094791e+02];% Coeficientes do polinômio da  bomba
        P.q0_min = 1000/3600; % Vazão mínima de referência
        P.q0_max = 4000/3600; % Vazão máxima de referência
        P.nestagios = 15;                   %        Número de estágios da bomba 
        P.f0 = 60;                          % [Hz]   Frequência de referência da curva da bomba
        %P.dfMAX = 0.05/P.Ts;                 % [Hz/s] Variação máxima da frequência 0.5Hz/s
        P.fMax = 60;                        % [Hz]   Frequência máxima admissível
        P.fMin = 30;                        % [Hz]   Frequência mínima admissível
        %P.dzMAX = 0.05/P.Ts;                 % [%/s]  Variação máxima admissível
        P.z_cMax = 1;                     % [%] Abertura máxima admissível
        P.z_cMin = 0;                       % [%] abertura mínima admissível
        
        f = 30:0.01:60; % Range de frequência da bomba
        H0_dt = P.Head(1)*P.q0_min^2 + P.Head(2)*P.q0_min + P.Head(3); % Curva de Head de referência máximo
        H0_ut = P.Head(1)*P.q0_max^2 + P.Head(2)*P.q0_max + P.Head(3); % Curva de Head de referência mínimo
        H_ut = H0_ut*(f./P.f0).^2; % Curva de Head
        H_dt = H0_dt*(f./P.f0).^2; % Curva de Head
        P.Hdown_min = H_dt(1); % Valor mínimo da curva de Downthrust
        P.HUp_min = H_ut(1); % Valor mínimo da curva de Uptrhrust
        %% Parâmetros de admensionalização do modelo
        load('Envelope_BCS_LEA.mat','Head_60Hz','Head_30Hz')
        P.DeltaP_max = ((max(Head_60Hz)-min(Head_30Hz))*P.rho(1)*P.g)/10^5;
        P.Delta_vazao = 1000*4/3600;
        P.L_poco = P.h(1);  
        %% Seleção das condições iniciais análogas aos dados experimentais coletados no LEA.

        P.z_c = 100;                                                   % [%]    Posição inicial da choke
        P.z_cRef = P.z_c;                                              % [%]    Referência de abertura da válvula, na condição inicial (setpoint)
        P.f_sp = 60;                                                   % [Hz]   Referência da frequência do BCS no início do experimento(setpoint)
        P.z_c_sp = 100;                                                % [%]    Referência de abertura da válvula (setpoint)
        P.Estado_inicial = [13.319668860033406,4.605391101170010,1.544104925333945*1000/3600];
        
        P.q   = P.Estado_inicial(3);             % [m³/s] Condição inicial da vazão média
        P.q_r = P.Estado_inicial(3);             % [m³/s] Condição inicial da vazão do reservatório
        P.q_p = P.Estado_inicial(3);             % [m³/s] Condição inicial da vazão da bomba
        P.q_c = P.Estado_inicial(3);             % [m³/s] Condição inicial da vazão na choke (produção)
        
        % Perfil de pressão na ordem: Pressão de reservatório, pressão de fundo,
        % pressão de intake, pressão de descarga, pressão na choke, pressão no
        % manifold e incremento de pressão da bomba (Delta P).
        P.p = [2.286875247955322e+05,1.897236035469732*10^5,1.117193885469732*10^5,6.509366267692657*10^5,4.605391101170010*10^5,6.200329065322876e+04,5.392172382222925*10^5]/(10^5);
           
           
        %% Parâmetros de admensionalização do modelo
        load('Envelope_BCS_LEA.mat','Head_60Hz','Head_30Hz')
        P.DeltaP_max = ((max(Head_60Hz)-min(Head_30Hz))*P.rho(1)*P.g)/10^5;
        P.Delta_vazao = 1000*4/3600;
        P.L_poco = P.h(1);       

end