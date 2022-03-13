function [Med,Saidas,D] = M_BCS(P,Entradas,D,n)
%%
Saidas = zeros(P.nsim,8);
Saidas(1,:) = [P.Estado_inicial(1)/P.L_poco P.Estado_inicial(2)/P.DeltaP_max P.Estado_inicial(3)/P.Delta_vazao P.q_c P.q_r P.p(2) P.p(3) P.p(7)];
Med = zeros(P.nsim,4);
Med(1,:) = [P.Estado_inicial(1) P.Estado_inicial(2) P.Estado_inicial(3) (P.q_c)*3.6];

switch n
    case 'malha aberta'
        %%
        h = waitbar(0,'Executando simulação por favor aguarde')
        for i=2:P.nsim
            waitbar(i/P.nsim,h)
            options = odeset('RelTol',1e-10,'Stats','off');
            Vetor_entradas = [Entradas.Frequencia(i) Entradas.Abertura_valvula_choke(i) Entradas.Pressao_reservatorio(i) Entradas.Pressao_manifold(i)];
            tic
            [~,y]=ode15s(@dinamica_BCS, [0 P.Ts], Saidas(i-1,1:3), options, Vetor_entradas,P);
            t(i)=toc;


            N_anular = y(end,1)*P.L_poco;
            P_wh = y(end,2)*P.DeltaP_max;
            q_media = y(end,3)*P.Delta_vazao;


            % Cálculo da pressão de intake.
            P_intake =   N_anular*P.rho(1)*P.g/10^5;
            if P_intake<0
                P_intake = 0;
            end 
            % Cálculo da pressão do fundo
            P_BHP = (P.rho(1)*P.g*(P.h(1)-P.h(2))+P.rho(1)*P.g*N_anular)/10^5;
            % Cálculo da vazão do reservatório
            q_reservatorio = ((1).*((Entradas.Pressao_reservatorio(i) - P_BHP).^0.5)*P.k_resv);
            % Vazão através da válvula choke.
            %q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            % Cálculo do incremento de pressão na bomba
            q_0a = q_media;%*(P.f0/Entradas.Frequencia(i));
            %Head corrigido
            H_0 = P.Head(1)*q_0a^2 +P.Head(2)*q_0a + P.Head(3);
            H(i) = H_0*(Entradas.Frequencia(i)/P.f0)^2;
            % Delta P sairá em Pascal
            DeltaP = H(i)*P.rho(2)*P.g/10^5;




            % As demais variáveis do sistema são armazenadas nesta matriz com a
            % seguinte ordem: 1 Vazão média na coluna; 2 vazão do reservatório; 3
            % Pressão no fundo do poço; 3 Pressão de intake e 4 Delta de pressão da
            % bomba
            %%%

            Saidas(i,:) = [y(end,1) y(end,2)  y(end,3) q_choke  q_reservatorio  P_BHP P_intake DeltaP];
             
        end
    close(h)     
    %%
     
    case 'malha aberta filtro'
        %%
%        h = waitbar(0,'Executando a simulação, por favor aguarde.');
        for i=2:P.nsim
%            waitbar(i/P.nsim,h);
            options = odeset('RelTol',1e-10,'Stats','off');
            Vetor_entradas = [Entradas.Frequencia(i) Entradas.Abertura_valvula_choke(i) Entradas.Pressao_reservatorio(i) Entradas.Pressao_manifold(i)];
            tic
            [~,y]=ode15s(@dinamica_BCS_filtro, [0 P.Ts], Saidas(i-1,1:3), options, Vetor_entradas,P);
            t(i)=toc;


            N_anular = y(end,1)*P.L_poco;
            P_wh = y(end,2)*P.DeltaP_max;
            q_media = y(end,3)*P.Delta_vazao;


            % Cálculo da pressão de intake.
            P_intake =   N_anular*P.rho(1)*P.g/10^5;
            if P_intake<0
                P_intake = 0;
            end 
            % Cálculo da pressão do fundo
            P_BHP = (P.rho(1)*P.g*(P.h(1)-P.h(2))+P.rho(1)*P.g*N_anular)/10^5;
            % Cálculo da vazão do reservatório
            q_reservatorio = ((1).*((Entradas.Pressao_reservatorio(i) - P_BHP).^0.5)*P.k_resv);
            % Vazão através da válvula choke.
            %q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            % Cálculo do incremento de pressão na bomba
            q_0a = q_media;%*(P.f0/Entradas.Frequencia(i));
            %Head corrigido
            H_0 = P.Head(1)*q_0a^2 +P.Head(2)*q_0a + P.Head(3);
            H(i) = H_0*(Entradas.Frequencia(i)/P.f0)^2;
            % Delta P sairá em Pascal
            DeltaP = H(i)*P.rho(2)*P.g/10^5;




            % As demais variáveis do sistema são armazenadas nesta matriz com a
            % seguinte ordem: 1 Vazão média na coluna; 2 vazão do reservatório; 3
            % Pressão no fundo do poço; 3 Pressão de intake e 4 Delta de pressão da
            % bomba
            %%%

            Saidas(i,:) = [y(end,1) y(end,2)  y(end,3) q_choke  q_reservatorio  P_BHP P_intake DeltaP];
%            close(h)   
        end
            
    
    case 'malha PID'
        %%
        D.PID.NIVEL.erro = zeros(P.nsim,1);
        D.PID.NIVEL.proporcional = zeros(P.nsim,1);
        D.PID.NIVEL.derivativo = zeros(P.nsim,1);
        D.PID.NIVEL.integrativo = zeros(P.nsim,1);
        
        D.PID.VAZAO.erro = zeros(P.nsim,1);
        D.PID.VAZAO.proporcional = zeros(P.nsim,1);
        D.PID.VAZAO.derivativo = zeros(P.nsim,1);
        D.PID.VAZAO.integrativo = zeros(P.nsim,1);
        
        D.simu.AF = zeros(P.nsim,1);
        D.simu.AA = zeros(P.nsim,1);
        D.simu.A = zeros(P.nsim,1);
        D.simu.A(1) = Entradas.Abertura_valvula_choke(1);
        D.simu.F = zeros(P.nsim,1);
        D.simu.F(1) = Entradas.Frequencia(1);
        
        D.PID.f_max = 60;
        D.PID.f_min = 30;
        D.PID.sf_max = 5;
        

        D.PID.a_max = 1;
        D.PID.a_min = 0.2;
        D.PID.sa_max = 0.2;
        
       
        D.PID.VAZAO.A = zeros(P.nsim,1);
        D.PID.VAZAO.A(1:2) = 90/100;
        D.PID.VAZAO.F = ones(P.nsim,1)*D.PID.f_max;
        D.PID.NIVEL.F = ones(P.nsim,1)*D.PID.f_max
        
        
        %h = waitbar(0,'Executando a simulação, por favor aguarde.');
            for i=2:P.nsim
         %       waitbar(i/P.nsim,h);
                D.PID.NIVEL.erro(i)= (D.PID.NIVEL.ref/P.L_poco - Saidas(i-1,1));
                D.PID.NIVEL.proporcional(i) = D.PID.NIVEL.erro(i);
                D.PID.NIVEL.derivativo(i) = ((D.PID.NIVEL.erro(i-1) - D.PID.NIVEL.erro(i))*P.Ts);
                D.PID.NIVEL.integrativo(i) = ((D.PID.NIVEL.erro(i-1)+(D.PID.NIVEL.erro(i)))/P.Ts);
                if(i<300)
                    D.PID.NIVEL.I = sum(D.PID.NIVEL.integrativo);
                elseif(i<600)
                    D.PID.NIVEL.I = sum(D.PID.NIVEL.integrativo(299:600));
                else
                    D.PID.NIVEL.I = sum(D.PID.NIVEL.integrativo(601:i-1));
                end

                SC_NF = -(D.PID.NIVEL.kp*D.PID.NIVEL.proporcional(i)+D.PID.NIVEL.kd*D.PID.NIVEL.derivativo(i)+ D.PID.NIVEL.ki*D.PID.NIVEL.I);
                D.simu.AF(i) = SC_NF;
                if(abs(SC_NF) > D.PID.sf_max)
                    D.PID.NIVEL.F(i) = D.PID.NIVEL.F(i-1) + D.PID.sf_max*(SC_NF/abs(SC_NF));
                else
                    D.PID.NIVEL.F(i) = D.PID.NIVEL.F(i-1) + SC_NF;
                end
                if(D.PID.NIVEL.F(i) > D.PID.f_max)
                    D.PID.NIVEL.F(i) = D.PID.f_max;
                end
                if(D.PID.NIVEL.F(i) < D.PID.f_min)
                    D.PID.NIVEL.F(i) = D.PID.f_min;
                end

                D.PID.VAZAO.erro(i)= (D.PID.VAZAO.ref*1000/3600 - Saidas(i-1,4));
                D.PID.VAZAO.proporcional(i) = D.PID.VAZAO.erro(i);
                D.PID.VAZAO.derivativo(i) = ((D.PID.VAZAO.erro(i-1) - D.PID.VAZAO.erro(i))*P.Ts);
                D.PID.VAZAO.integrativo(i) = ((D.PID.VAZAO.erro(i-1)+(D.PID.VAZAO.erro(i)))/P.Ts);
                if(i<300)
                    D.PID.VAZAO.I = sum(D.PID.VAZAO.integrativo);
                elseif(i<600)
                    D.PID.VAZAO.I = sum(D.PID.VAZAO.integrativo(299:600));
                else
                    D.PID.VAZAO.I = sum(D.PID.VAZAO.integrativo(601:i-1));
                end
                %if(i == m_override)
                %    D.PID.VAZAO.F_D.PID(i-1) = D.PID.NIVEL.F_D.PID(i-1);
                %end
                %Seletor de override para vazão em razão de frequencia


                if(D.PID.override.modo == 0)
                SC_VA = (D.PID.VAZAO.kp*D.PID.VAZAO.proporcional(i)+D.PID.VAZAO.kd*D.PID.VAZAO.derivativo(i)+ D.PID.VAZAO.ki*D.PID.VAZAO.I);
                D.simu.AA(i) = SC_VA;
                if(abs(SC_VA) > D.PID.sa_max)
                    D.PID.VAZAO.A(i) = D.PID.VAZAO.A(i-1) + D.PID.sa_max*(SC_VA/abs(SC_VA));
                else
                    D.PID.VAZAO.A(i) = D.PID.VAZAO.A(i-1) + SC_VA;
                end

                if(D.PID.VAZAO.A(i) > D.PID.a_max)
                    D.PID.VAZAO.A(i) = D.PID.a_max;
                end
                if(D.PID.VAZAO.A(i) < D.PID.a_min)
                    D.PID.VAZAO.A(i) = D.PID.a_min;
                end
                Entradas.Frequencia(i) = D.PID.NIVEL.F(i);
                
                elseif(D.PID.override.modo == 1)
                    %em caso de override 

                    D.PID.VAZAO.A(i) = 1;
                    SC_VA = (D.PID.VAZAO.kpf*D.PID.VAZAO.proporcional(i)+D.PID.VAZAO.kdf*D.PID.VAZAO.derivativo(i)+ D.PID.VAZAO.kif*D.PID.VAZAO.I);
                    D.simu.AA(i) = SC_VA;
                    if(abs(SC_VA) > D.PID.sf_max)
                        D.PID.VAZAO.F(i) = D.PID.VAZAO.F(i-1) + D.PID.sf_max*(SC_VA/abs(SC_VA));
                    else
                        D.PID.VAZAO.F(i) = D.PID.VAZAO.F(i-1) + SC_VA;
                    end

                    if(D.PID.VAZAO.F(i) > D.PID.f_max)
                        D.PID.VAZAO.F(i) = D.PID.f_max;
                    end
                    if(D.PID.VAZAO.F(i) < D.PID.f_min)
                        D.PID.VAZAO.F(i) = D.PID.f_min;
                    end
                    if(D.PID.NIVEL.F(i) <= D.PID.VAZAO.F(i))
                        Entradas.Frequencia(i) = D.PID.NIVEL.F(i);
                    else
                        Entradas.Frequencia(i) = D.PID.VAZAO.F(i);
                    end
                end
                Entradas.Abertura_valvula_choke(i) = D.PID.VAZAO.A(i);
                Entradas.Frequencia(i) = D.PID.NIVEL.F(i);
                D.simu.A(i) = Entradas.Abertura_valvula_choke(i);
                D.simu.F(i) = Entradas.Frequencia(i);
        %%% Fim do controlador

        options = odeset('RelTol',1e-10,'Stats','off');
        Vetor_entradas = [(Entradas.Frequencia(i)- (Entradas.Frequencia(i-1)*0.1))  Entradas.Abertura_valvula_choke(i) Entradas.Pressao_reservatorio(i) Entradas.Pressao_manifold(i)];
        tic
        [~,y]=ode15s(@dinamica_BCS, [0 P.Ts], Saidas(i-1,1:3), options, Vetor_entradas,P);
        t(i)=toc;
        
         ruido = randn(3,1)*0.001;
         ruido(1) = ruido(1)*0.1;
         y(end,1) = y(end,1) + ruido(1);
         y(end,2) = y(end,2) + ruido(2);
         y(end,3) = y(end,3) + ruido(3);
        
        N_anular = y(end,1)*P.L_poco;
        P_wh = y(end,2)*P.DeltaP_max;
        q_media = y(end,3)*P.Delta_vazao;


        % Cálculo da pressão de intake.
        P_intake =   N_anular*P.rho(1)*P.g/10^5;
        if P_intake<0
            P_intake = 0;
        end 
        % Cálculo da pressão do fundo
        P_BHP = (P.rho(1)*P.g*(P.h(1)-P.h(2))+P.rho(1)*P.g*N_anular)/10^5;
        % Cálculo da vazão do reservatório
        q_reservatorio = ((1).*((Entradas.Pressao_reservatorio(i) - P_BHP).^0.5)*P.k_resv);
        % Vazão através da válvula choke.
        %q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
        q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
        % Cálculo do incremento de pressão na bomba
        q_0a = q_media;%*(P.f0/Entradas.Frequencia(i));
        %Head corrigido
        H_0 = P.Head(1)*q_0a^2 +P.Head(2)*q_0a + P.Head(3);
        H(i) = H_0*(Entradas.Frequencia(i)/P.f0)^2;
        % Delta P sairá em Pascal
        DeltaP = H(i)*P.rho(2)*P.g/10^5;

        
 
        
        % As demais variáveis do sistema são armazenadas nesta matriz com a
        % seguinte ordem: 1 Vazão média na coluna; 2 vazão do reservatório; 3
        % Pressão no fundo do poço; 3 Pressão de intake e 4 Delta de pressão da
        % bomba
        %%%

        Saidas(i,:) = [y(end,1) y(end,2)  y(end,3) q_choke  q_reservatorio  P_BHP P_intake DeltaP];
            end
        %close(h)
    %%
    case 'malha aberta com ruido'
       cont = 0;
       h = waitbar(0,'Carregando simulação');
       D.entradas.frequencia = zeros(length(Entradas.Frequencia),1);
       D.entradas.frequencia_desvio = D.entradas.frequencia;
       D.entradas.abertura = D.entradas.frequencia;
       D.entradas.abertura_desvio = D.entradas.frequencia;
       
       D.entradas.frequencia(1) = Entradas.Frequencia(1);
       D.entradas.frequencia_desvio(1) = D.entradas.frequencia(1)- D.Estacionario.frequencia;
       D.entradas.abertura(1) = Entradas.Abertura_valvula_choke(1);
       D.entradas.abertura_desvio(1) = D.entradas.abertura(1) - D.Estacionario.abertura_valvula_choke;
       for i=2:P.nsim
            cont = cont +1;
            fprintf(['Contador:',num2str(cont),'\n'])
            waitbar(i/P.nsim,h);
            % waitbar(i/P.nsim,h);
            options = odeset('RelTol',1e-10,'Stats','off');
            erro = 0;
            
            D.entradas.frequencia(i) = Entradas.Frequencia(i)+randn(1)*10^-5;
            D.entradas.frequencia_desvio(i) = D.entradas.frequencia(i) - D.Estacionario.frequencia;
            D.entradas.abertura(i) = Entradas.Abertura_valvula_choke(i)+randn(1)*10^-5;
            D.entradas.abertura_desvio(i) = D.entradas.abertura(i) - D.Estacionario.abertura_valvula_choke;
            
            
            Vetor_entradas = [Entradas.Frequencia(i)+D.Fa(i,1) Entradas.Abertura_valvula_choke(i)+D.Fa(i,2) (Entradas.Pressao_reservatorio(i) - erro) Entradas.Pressao_manifold(i)];
            tic
            [~,y]=ode15s(@dinamica_BCS, [0 P.Ts], Saidas(i-1,1:3), options, Vetor_entradas,P);
            t(i)=toc;
            
            N_anular = y(end,1)*P.L_poco;
            P_wh = y(end,2)*P.DeltaP_max;
            q_media = y(end,3)*P.Delta_vazao;
            
            % Cálculo da pressão de intake.
            P_intake =   N_anular*P.rho(1)*P.g/10^5;
            if P_intake<0
                P_intake = 0;
            end 
            % Cálculo da pressão do fundo
            P_BHP = (P.rho(1)*P.g*(P.h(1)-P.h(2))+P.rho(1)*P.g*N_anular)/10^5;
            % Cálculo da vazão do reservatório
            q_reservatorio = ((1).*((Entradas.Pressao_reservatorio(i) - P_BHP).^0.5)*P.k_resv);
            % Vazão através da válvula choke.
            %q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            % Cálculo do incremento de pressão na bomba
            q_0a = q_media;%*(P.f0/Entradas.Frequencia(i));
            %Head corrigido
            H_0 = P.Head(1)*q_0a^2 +P.Head(2)*q_0a + P.Head(3);
            H(i) = H_0*(Entradas.Frequencia(i)/P.f0)^2;
            % Delta P sairá em Pascal
            DeltaP = H(i)*P.rho(2)*P.g/10^5;

            % As demais variáveis do sistema são armazenadas nesta matriz com a
            % seguinte ordem: 1 Vazão média na coluna; 2 vazão do reservatório; 3
            % Pressão no fundo do poço; 3 Pressão de intake e 4 Delta de pressão da
            % bomba
            %%%
            
            ruido = randn(4,1)*D.pertu.medicao_var;
            %ruido(4) = ruido(4)*10;
            ruido(4) = ruido(1)/10;
            %D.Saida.nivel
            Saidas(i,:) = [(y(end,1)) y(end,2)  (y(end,3)) q_choke  q_reservatorio  P_BHP P_intake DeltaP];
            %close(h)
            Med(i,:) = [(y(end,1))*P.L_poco+ruido(1)+D.Fsy(i,1) (y(end,2))*P.DeltaP_max+ruido(2)+D.Fsy(i,2) (y(end,3))*P.Delta_vazao+ruido(3) (q_choke)+ruido(4)+D.Fsy(i,3)];
       end
        
        
       close(h)
    case 'linearização'
        for i=2:P.nsim
%            waitbar(i/P.nsim,h);
            
            

            options = odeset('RelTol',1e-10,'Stats','off');
            Vetor_entradas = [D.equi.F D.equi.A D.equi.Pr D.equi.Pm];
            tic
            [~,y]=ode15s(@dinamica_BCS, [0 P.Ts], Saidas(i-1,1:3), options, Vetor_entradas,P);
            t(i)=toc;


            N_anular = y(end,1)*P.L_poco;
            P_wh = y(end,2)*P.DeltaP_max;
            q_media = y(end,3)*P.Delta_vazao;


            % Cálculo da pressão de intake.
            P_intake =   N_anular*P.rho(1)*P.g/10^5;
            if P_intake<0
                P_intake = 0;
            end 
            % Cálculo da pressão do fundo
            P_BHP = (P.rho(1)*P.g*(P.h(1)-P.h(2))+P.rho(1)*P.g*N_anular)/10^5;
            % Cálculo da vazão do reservatório
            q_reservatorio = ((1).*((Entradas.Pressao_reservatorio(i) - P_BHP).^0.5)*P.k_resv);
            % Vazão através da válvula choke.
            %q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            % Cálculo do incremento de pressão na bomba
            q_0a = q_media;%*(P.f0/Entradas.Frequencia(i));
            %Head corrigido
            H_0 = P.Head(1)*q_0a^2 +P.Head(2)*q_0a + P.Head(3);
            H(i) = H_0*(Entradas.Frequencia(i)/P.f0)^2;
            % Delta P sairá em Pascal
            DeltaP = H(i)*P.rho(2)*P.g/10^5;




            % As demais variáveis do sistema são armazenadas nesta matriz com a
            % seguinte ordem: 1 Vazão média na coluna; 2 vazão do reservatório; 3
            % Pressão no fundo do poço; 3 Pressão de intake e 4 Delta de pressão da
            % bomba
            %%%

            Saidas(i,:) = [y(end,1) y(end,2)  y(end,3) q_choke  q_reservatorio  P_BHP P_intake DeltaP];
        end
    case 'malha aberta falha'
        ruido = randn(3,1)*0;
        h = waitbar(0,'Executando simulação por favor aguarde')
        fprintf('o nivel da falha é [%2.1f %2.1f]\n',D.falha(1,1),D.falha(1,2))
        for i=2:P.nsim
            waitbar(i/P.nsim,h)
            options = odeset('RelTol',1e-10,'Stats','off');
            Vetor_entradas = [Entradas.Frequencia(i) Entradas.Abertura_valvula_choke(i) Entradas.Pressao_reservatorio(i) Entradas.Pressao_manifold(i)];
            tic
            [~,y]=ode15s(@dinamica_BCS_falha, [0 P.Ts], Saidas(i-1,1:3), options, Vetor_entradas,P,D.falha(i,:));
            t(i)=toc;
            if ((D.falha(i-1,1) ~= D.falha(i,1)) ||  (D.falha(i-1,2) ~= D.falha(i,2)))
                fprintf('o nivel da falha é [%2.1f %2.1f] no momento %i\n',D.falha(i,1),D.falha(i,2),i*P.Ts)
            end
            
            N_anular = y(end,1)*P.L_poco;
            P_wh = y(end,2)*P.DeltaP_max;
            q_media = y(end,3)*P.Delta_vazao;


            % Cálculo da pressão de intake.
            P_intake =   N_anular*P.rho(1)*P.g/10^5;
            if P_intake<0
                P_intake = 0;
            end 
            % Cálculo da pressão do fundo
            P_BHP = (P.rho(1)*P.g*(P.h(1)-P.h(2))+P.rho(1)*P.g*N_anular)/10^5;
            % Cálculo da vazão do reservatório
            q_reservatorio = ((1).*((Entradas.Pressao_reservatorio(i) - P_BHP).^0.5)*P.k_resv)-0.2*D.falha(1);
            % Vazão através da válvula choke.
            %q_choke = ((Entradas.Abertura_valvula_choke(i)).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            q_choke = ((D.falha(2)+Entradas.Abertura_valvula_choke(i)*(1-D.falha(2))).*((P_wh-Entradas.Pressao_manifold(i)).^0.5)*P.k_choke);
            % Cálculo do incremento de pressão na bomba
            q_0a = q_media;%*(P.f0/Entradas.Frequencia(i));
            %Head corrigido
            H_0 = P.Head(1)*q_0a^2 +P.Head(2)*q_0a + P.Head(3);
            H(i) = H_0*(Entradas.Frequencia(i)/P.f0)^2;
            % Delta P sairá em Pascal
            DeltaP = H(i)*P.rho(2)*P.g/10^5;




            % As demais variáveis do sistema são armazenadas nesta matriz com a
            % seguinte ordem: 1 Vazão média na coluna; 2 vazão do reservatório; 3
            % Pressão no fundo do poço; 3 Pressão de intake e 4 Delta de pressão da
            % bomba
            %%%
            Saidas(i,:) = [y(end,1) y(end,2)  y(end,3) q_choke  q_reservatorio  P_BHP P_intake DeltaP];
            Med(i,:) = [N_anular+(D.Estacionario.nivel*0.001)*ruido(1) P_wh+(D.Estacionario.pressao_choke*0.01)*ruido(2) q_choke+(D.Estacionario.vazao_choke*0.01)*ruido(3) q_choke*3.6+(D.Estacionario.vazao_choke*0.01)*ruido(3)];
        end
        close(h)
    otherwise
            error('Defina o case corretamente.')
    end
end