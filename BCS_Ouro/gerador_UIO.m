function [M] = gerador_UIO(ssx,E,dados)
    function [C,Ceq] = NON(x,ss)
        C = eig(ss.A - x*ss.C);
        Ceq = [];
    end
    % Condição de relação dos disturbios com as saidas
    if(rank(ssx.C*E) ~= rank(E))
        error("rank de CE(%i) é diferente que rank de E(%i)\n",rank(ssx.C*E),rank(E))
    end
    M.H = E*pinv(ssx.C*E); % H = E*(CE)+
    M.T = eye(length(M.H*ssx.C))- M.H*ssx.C; % T = I - H*C
    M.TB = M.T*ssx.B;% T*B
    A1 = M.T*ssx.A; % A*T
    f=@(K)sum(eig(A1 - K*ssx.C)) % autovalores de F
    t_K = [size(A1,2) size(ssx.C,1)]; 
    [M.K1,~,~,~] = fmincon(f,ones(t_K),[],[],[],[],ones(t_K),dados.K*ones(t_K),@(x)NON(x,ssx));
    M.F = A1 - M.K1*ssx.C; % F = TA - K1*C
    fprintf('autovalores:%d\n',eig(M.F))
    M.K2 = M.F*M.H; % K2 = F*H
    M.K = M.K1+M.K2; % K = K1+K2
    % Observador
    M.sszx = ss(M.F,[M.TB M.K],eye(length(M.F)),[zeros(size(M.F,1),size(M.TB,2)) M.H]);
    M.sszx.OutputName = {'x'};
    M.sszy = ss(M.F,[M.TB M.K],ssx.C,[zeros(size(ssx.C,1),size(M.TB,2)) ssx.C*M.H]);
    
    

end