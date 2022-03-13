function [P,dados] = estacionario(entrada,P,dados)

    Entradas.Frequencia = entrada(1)*ones(P.nsim,1);
    Entradas.Abertura_valvula_choke = entrada(2)*ones(P.nsim,1);
    Entradas.Pressao_reservatorio = P.p(1)*ones(length(Entradas.Frequencia),1);
    Entradas.Pressao_manifold = P.p(6)*ones(length(Entradas.Frequencia),1);
    
    fprintf('    Levando o sistema para o estado estacionário desejado.\n\n')
    n = 'malha aberta';
    n_a = 'LEA';
    [~,Saidas,dados] = M_BCS(P,Entradas,dados,n);


    %Entradas.Pressao_reservatorio = P.p(1);
    %Entradas.Pressao_manifold = P.p(6);
    opt=optimset('TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',10000,'Display','off');
    [xopt,fval,flag,iter]=fsolve(@(x)dinamica_BCS(1,x,[Entradas.Frequencia(1) Entradas.Abertura_valvula_choke(1) P.p(1) P.p(6)],P),Saidas(end,1:3),opt);

    if flag >0
        fprintf('    O ponto encontrado é um estado estacionário.\n\n')
    elseif flag==0 
        fprintf('    O ponto especficado pode não ser um estado estacionário. Aumente o tempo de simulação da pré-solução para confirmar.\n\n')
    else
        error('Não foi possível encontrar um ponto estacionário. Reveja os parâmetros do sistema.')
    end
        dados.Estacionario.nivel = xopt(1).*P.L_poco; % Metro
        dados.Estacionario.pressao_choke = xopt(2).*P.DeltaP_max; % Bar
        dados.Estacionario.vazao_media =  xopt(3).*3.6*P.Delta_vazao; % m³/h
        dados.Estacionario.pressao_intake = dados.Estacionario.nivel*P.rho(1)*P.g/10^5;
        dados.Estacionario.pressao_fundo = (dados.Estacionario.pressao_intake*10^5 + P.rho(1)*P.g*(P.h(1)-P.h(2)))/10^5;
        dados.Estacionario.vazao_reservatorio = ((1).*((P.p(1) - dados.Estacionario.pressao_fundo).^0.5)*P.k_resv)*3.6;
        dados.Estacionario.vazao_choke = ((Entradas.Abertura_valvula_choke(1)).*((dados.Estacionario.pressao_choke-P.p(6)).^(0.5))*P.k_choke)*3.6;
        dados.Estacionario.DeltaP = dinamica_dp(dados.Estacionario.vazao_choke/3600,Entradas.Frequencia(1),P,n_a)/10^5;
        dados.Estacionario.Head = (dados.Estacionario.DeltaP*10^5)/(P.g*P.rho(1));
        dados.Estacionario.pressao_descarga_bomba = dados.Estacionario.pressao_intake + dados.Estacionario.DeltaP;
        dados.Estacionario.frequencia = Entradas.Frequencia(1);
        dados.Estacionario.abertura_valvula_choke = Entradas.Abertura_valvula_choke(1);
        dados.Estacionario.pressao_manifold = P.p(6);
        dados.Estacionario.pressao_reservatorio = P.p(1);
        
   
        
        % Perfil de pressão na ordem: Pressão de reservatório, pressão de fundo,
        % pressão de intake, pressão de descarga, pressão na choke, pressão no
        % manifold e incremento de pressão da bomba (Delta P).
        
        P.p = [P.p(1) dados.Estacionario.pressao_fundo dados.Estacionario.pressao_intake ...
               dados.Estacionario.pressao_descarga_bomba dados.Estacionario.pressao_choke...
               P.p(6) dados.Estacionario.DeltaP ];
        P.Estado_inicial = xopt;
        P.q   = xopt(3)*P.Delta_vazao;             % [L/s] Condição inicial da vazão média
        P.q_r = xopt(3)*P.Delta_vazao;             % [L/s] Condição inicial da vazão do reservatório
        P.q_p = xopt(3)*P.Delta_vazao;             % [L/s] Condição inicial da vazão da bomba
        P.q_c = xopt(3)*P.Delta_vazao;             % [L/s] Condição inicial da vazão na choke (produção)
        
        P.Estado_inicial = [xopt(1)*P.L_poco xopt(2)*P.DeltaP_max xopt(3)*P.Delta_vazao];
end