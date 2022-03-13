function [M] = gerador_UIO_reduced(sys,E,i,R,Q,dados)
    M.D = sys.D;
      
    function [C,Ceq] = NON(x,A,C)
            C = eig(A - x*C);
            Ceq = [];
    end

    if(rank(sys.C*E(:,i)) == rank(E(:,i)))
        fprintf('Condição 1 aceita\n')
    else
        error('problema na condição 1')
    end
    
    M.ordem = rank(E(:,i));
    [U1,S1,V1] = svd(E(:,i));
    
    M.T = transpose(U1);
    TAT = M.T*sys.A*inv(M.T);
    A11 = TAT(1,1);
    A12 = TAT(1,2:end);
    A21 = TAT(2:end,1);
    A22 = TAT(2:end,2:end);
    
    TB = M.T*sys.B;
    M.CT = sys.C*inv(M.T);
    TE1 = M.T*E(:,i);
    TE2 = M.T*[E(:,1:i-1),E(:,i+1:end)];
    [U2,S2,V2] = svd(M.CT(:,1));
    Sigma2 = S2(1);
    M.T1 = transpose(U2);
    UC = M.T1*M.CT(:,2:end);
    M.C21 = UC(1,:);
    C22 = UC(2:end,:);
    B1 = TB(1,:);
    B2 = TB(2:end,:);
    A22_barra = A22-A21*pinv(Sigma2*V2')*M.C21;
    M.ssz = ss(A22_barra,[B2 A21*(Sigma2*V2') zeros(size(A22,1),rank(sys.C)-rank(E(:,i)))],eye(length(A22)),[0])
    
    M.ssy = ss(A22_barra,[B2 A21*(Sigma2*V2') zeros(size(A22,1),rank(sys.C)-rank(E(:,i)))],C22,[0])

    M.C22 = C22;
    
    M.z1 =@(z2,y1)pinv(Sigma2*V2')*(y1' - M.C21*z2')';
    M.f1 =@(z2,y1,u) pinv(TE1(1))*(M.z1(z2(2,:),y1(2))' - A11*M.z1(z2(1,:),y1(1))' -A12*z2(1,:)' -B1*u')';
    
    f=@(K,A,C)sum(eig(A - K*C));
    t_K = [size(M.ssy.A,1) size(M.ssy.C,1)]; 
    [M.K,M.AV,~,~] = fmincon(@(x)f(x,M.ssy.A,M.ssy.C),ones(t_K),[],[],[],[],zeros(t_K),dados.K*ones(t_K),@(x)NON(x,M.ssy.A,M.ssy.C));
    if(min(M.AV) < 0 )
        error('Observador instaveis')
    end
    M.ssy.B(:,size([B2 A21*(Sigma2*V2')],2)+1:end) = M.K;
    M.ssz.B(:,size([B2 A21*(Sigma2*V2')],2)+1:end) = M.K;
    M.s_obs = size([B2 A21*(Sigma2*V2')]);
    
    t_B = size(sys.B,2);
    t_y1 = rank(E(:,i));
    t_y2 = size(sys.C,1) - t_y1; 
    name = cell(1);
    for k = 1:t_B
        name{k} = {['u(',(num2str(k)),')']};
        M.ssy.InputName(k) = {['u(',(num2str(k)),')']};
    end
    for k = t_B+1:t_B+t_y1
        name{k} = {['y1(',(num2str(k-t_B)),')']};
        M.ssy.InputName(k) = {['y1(',(num2str(k-t_B)),')']};
    end
    for k = t_B+t_y1+1:t_B+t_y1+t_y2
        name{k} =  {['r(',(num2str(k-(t_B+t_y1))),')']};
        M.ssy.InputName(k) = {['r(',(num2str(k-(t_B+t_y1))),')']};
    end
    
    M.t_y1 = t_y1;
    M.t_y2 = t_y2;
    M.ssz.OutputName = {'z2'};
    %M.ssy.StateName = {'z'};
    %M.ssy.InputName =  {'u(1)','u(2)','y1','r(1)','r(2)'};
    M.ssy.OutputName = {'y2_h'};
    M.S1 = sumblk('r = y2 - y2_h',size(sys.C,1) - t_y1);
    
    %M.ss = connect(M.ssy,M.S1,...
     %   {'u(1)','u(2)','y1','y2(1)','y2(2)'},{'r(1)','r(2)','y2_h'});
    M.ss = connect(M.ssy,M.S1,{'u','y1','y2'},{'r','y2_h'});
        
    M.ss.StateName = {'z2(1)','z2(2)'};
    
    M.ssx = ss(M.ssy.A,[M.ssy.B(:,1:3) M.ssy.B(:,1:3)],M.ssy.C,[0])
    
    [M.KF,~,~] = kalman(M.ssx,Q,R);
    
    if(rank(obsv(M.ssy))== rank(M.ssy.A))
        fprintf('Condição 2 aceita\n')
    else
        error('problema na condição 2')
    end
end