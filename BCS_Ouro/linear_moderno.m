
function [P,dados] = linear_moderno(entrada,P,dados)
        %entrada = [45 0.95];
        [P,dados] = estacionario(entrada,P,dados);
        %%
        syms  Na Qm Pc f Z fref zref H Qc a1 a2


%%

        Pm = dados.Estacionario.pressao_manifold;
        Pr = dados.Estacionario.pressao_reservatorio;
        %% Equações algébricas
        h_c =  P.h(3);                            % [m]      Profundidade do choke
        h_p =  P.h(2);                            % [m]      Profundidade do conjunto BCS
        h_r =  P.h(1);                            % [m]      Profundidade do reservatório
%         L_poco = P.L_poco;
%         DeltaP_max = P.DeltaP_max;
%         Delta_vazao = P.Delta_vazao;
        Pf = (P.rho(1)*P.g*(P.h(1)-P.h(2)+Na))/10^5; % Pressão no fundo do poço
        Qr = 1.*P.k_resv(1).*(Pr - Pf).^(0.5); % Vazão do reservatório
        Qc = ((Z+a2).*((Pc-Pm).^0.5)*P.k_choke(1)); % Vazão na choke
        Q0 = Qm*(P.f0/f); % Vazão de referência do polinômio
        H0 = P.Head(1)*(Q0)^2 + P.Head(2)*(Q0) + P.Head(3); % Altura manométrica de referência
        H = H0*((f/P.f0)^2); % Correção do Head com leis de afinidade
        Dp = P.rho(1)*P.g*H; % Diferencial de pressão da bomba
        Re1 = 4*P.rho(1)*(Qm)/(pi*P.mu(1))/1000;
        Re2 = 4*P.rho(2)*(Qm)/(pi*P.mu(2))/1000;
        F_fric_1 =64/Re1;
        F_fric_2 =64/Re2;
        DF1 = F_fric_1*P.rho(1)*(((Qm)/1000)^2)/(2*pi*P.r(1)^3);
        DF2 = F_fric_2*P.rho(2)*(((Qm)/1000)^2)/(2*pi*P.r(2)^3);
        DF = DF1 + DF2;
        P.dtf = P.Ts;
        P.dtz = P.Ts;
        %% Construção das equações diferenciais
        % Equações diferenciais idênticas às existentes no arquivo "dinamica_bcs.m"
        f1 = (Qr + a1 - Qm)/(P.Aanular*1000); % Equação diferencial do nível
        f2 = ((P.B(2)/P.V(2))*(Qm-Qc)/10^3); % Equação diferencial da Pressão na choke
        f3 = (1000*(P.A_bar/(P.rho(1)*P.l_bar))*(Pf*10^5 - Pc*10^5 + Dp - DF -(P.rho(1)*P.g*(P.h(1)-P.h(2)))-(P.rho(2)*P.g*(P.h(2)-P.h(3))))); % Equação diferencial da vazão
        %f4 = (fref - f)/P.dtf; % Equação diferencial da resposta do inversor de freequêcia
        %f5 = (zref - Z)/P.dtz; % Equação diferencial da resposta da válvula choke

         %funcao=[f1;f2;f3;f4;f5];
         %funcao=[f1/P.L_poco;f2/P.DeltaP_max;f3/P.Delta_vazao];
         funcao=[f1;f2;f3];
         dados.F = funcao;
         
             %% Definição dos estados modelo 3x2
         
         Estados_3x2 = [Na Pc Qm];
         Entradas_3x2 = [f Z];
         Saidas_3x2 = [Na Pc Qc];
         falhas = [a1 a2];
         
         % Linearização
        A = jacobian(funcao,Estados_3x2);
        B = jacobian(funcao,Entradas_3x2);
        C = jacobian(Saidas_3x2,Estados_3x2);
        D = jacobian(Saidas_3x2,Entradas_3x2);
        E = jacobian(funcao,falhas);
        
        
        % Definição do ponto de linearização

        Z  = dados.Estacionario.abertura_valvula_choke;
        f  = dados.Estacionario.frequencia;
        Na = dados.Estacionario.nivel;
        Pc = dados.Estacionario.pressao_choke;
        Qm = dados.Estacionario.vazao_media*1000/3600;
        a1 = 0;
        a2 = 0;
        
        
        %dados.Estacionario = dados.Estacionario;
        dados.ss.simbolico.A = A;
        dados.ss.simbolico.B = B;
        dados.ss.simbolico.C = C;
        dados.ss.simbolico.D = D;
        dados.ss.simbolico.E = E;
        
        dados.ss.A = eval(A);
        dados.ss.B = eval(B);
        dados.ss.C = eval(C);
        dados.ss.D = eval(D);
        dados.ss.E = eval(E);
        dados.ss.modelo_linear = ss(eval(A),eval(B),eval(C),eval(D));
        
        %dados.ss.m_linear_falha = ss(eval(A),[eval(B) eval(E)],eval(C),eval(D));
        %% Normalização
        %{
        % x_dot = Gamma*A*Gamma^-1*
        Est = [dados.Estacionario.nivel dados.Estacionario.pressao_choke dados.Estacionario.vazao_media];
        dados.v_normal.A = [P.L_poco-Est(1),P.DeltaP_max-Est(2),P.Delta_vazao*3.6-Est(3)];
        dados.v_normal.C = [P.L_poco-Est(1),P.DeltaP_max-Est(2),P.Delta_vazao*3.6-dados.Estacionario.vazao_choke];
        G = inv(diag(dados.v_normal.A));
        Gy = inv(diag(dados.v_normal.C));
        dados.G = G;
        dados.Gy = Gy;
        %dados_linear.ss.ss.A = [dados_linear.ss.ss.A(1,:)/dados_linear.v_normal.A(1);dados_linear.ss.ss.A(2,:)/dados_linear.v_normal.A(2);dados_linear.ss.ss.A(3,:)/dados_linear.v_normal.A(3)];
        %dados_linear.ss.ss.B = [dados_linear.ss.ss.B(1,:)/dados_linear.v_normal.A(1);dados_linear.ss.ss.B(2,:)/dados_linear.v_normal.A(2);dados_linear.ss.ss.B(3,:)/dados_linear.v_normal.A(3)];
        %dados_linear.ss.ss.C = [dados_linear.ss.ss.C(1,:)/dados_linear.v_normal.A(1);dados_linear.ss.ss.C(2,:)/dados_linear.v_normal.A(2);dados_linear.ss.ss.C(3,:)/dados_linear.v_normal.A(3)];
        %dados_linear.ss.ss.D = [dados_linear.ss.ss.D(1,:)/dados_linear.v_normal.A(1);dados_linear.ss.ss.D(2,:)/dados_linear.v_normal.A(2);dados_linear.ss.ss.D(3,:)/dados_linear.v_normal.A(3)];
        %dados_linear.ss.E = [dados_linear.ss.E(1,:)/dados_linear.v_normal.A(1);dados_linear.ss.E(2,:)/dados_linear.v_normal.A(2);dados_linear.ss.E(3,:)/dados_linear.v_normal.A(3)];
        
        dados.ss.norm.A = G*dados.ss.A*inv(G);
        dados.ss.norm.B = G*dados.ss.B;
        dados.ss.norm.E = G*dados.ss.E;
        dados.ss.norm.C = Gy*dados.ss.C*inv(G);
        dados.ss.norm.D = Gy*dados.ss.D;
        dados.ss.norm.ss = ss(dados.ss.A,dados.ss.B,dados.ss.C,dados.ss.D);
        %eye(length(A)*0)
        %dados_linear.ss.dss = idss(eval(A),eval(B),eval(C),eval(D),[0],ones(length(A),1)*0,P.Ts);
        dados.ss.dss1 = c2d(ss(eval(A),eval(B),eval(C),eval([D])),P.Ts);
        dados.ss.dss2 = c2d(ss(eval(A),eval(E),eye(length(A)),[0]),P.Ts);

        dados.tf = dados.Modelo_linear.Unidades_engenharia.Expandido_3x2.Funcao_transferencia.FT;
        %}
         fprintf('\n Terminou \n')
end