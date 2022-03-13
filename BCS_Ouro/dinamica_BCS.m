%% Rotina que define o sistema de Equa��es Diferenciais Ordinarias do sistema BCS

function [dydt]=dinamica_BCS(t,y,Entradas,P)
% =========================================================================
% Define as entradas 
% =========================================================================
% Frequ�ncia
% Abertura da choke
% Press�o de reservat�rio
% Press�o de manifold

% =========================================================================
% Define os estados 
% =========================================================================
h0_anular = y(1)*P.L_poco; % Metros
P_wh = y(2)*P.DeltaP_max; % Bar
q_media = y(3)*P.Delta_vazao; % m^3/h
% fref = y(4);
% zref = y(5);
% =========================================================================
% Define a din�mica da v�lvula choP.ke, Fric��o e Deltap da bomba, Pressao de
% fundo e Vaz�o de reserv�torio
% =========================================================================



% C�lculo da perda de carga por fric��o.
Re = (4.*P.rho.*q_media)./(0.219*pi.*P.mu);

% =========================================================================
% Define o c�lculo da c�lculo da fric��o
% =========================================================================
if Re(1)<4000
fric = 64./Re;
else
fric = 0.36*Re.^(-0.25);
end

% =========================================================================
% Define o c�lculo da fri��o real
% =========================================================================
% A fric��o sai com unidade de kg/(ms�) equivalente a Bar
F = (fric.*P.rho.*q_media.^2)./(2.*pi.*P.r.^3);

H_aguabep = 4.330800000000000e+02; %ft
Q_aguabep = 4.401900000000000e+02; %bpd
yy = -112.1374+6.6504*log(H_aguabep)+12.8429*log(Q_aguabep);
Q = exp((39.5276+26.5605*log(P.mu(1)*1000)-yy)/51.6565); %Pa.s to Cp;
Cq = (1.0-4.0327e-03*Q-1.7240e-04*Q^2);
CH = 1.0-4.4723e-03*Q -4.1800e-05*Q^2; % 80%
q_0a = q_media/Cq*(P.f0/Entradas(1));
% q_0a = q_media;%*(P.f0/Entradas(1));
%Head corrigido
H_0 = P.Head(1)*q_0a^2 +P.Head(2)*q_0a + P.Head(3);
H = CH.*H_0*(Entradas(1)/P.f0)^2;
% Delta P sair� em Pascal
dp = H*P.rho(2)*P.g;
Qc = ((Entradas(2)).*((P_wh-Entradas(4)).^0.5)*P.k_choke(1));
P_BHP = (P.rho(1)*P.g*(P.h(1)-P.h(2))+P.rho(1)*P.g*h0_anular)./10.^5;
% Vaz�o no reservat�rio
if(Entradas(3)- P_BHP <= 0)
    fprintf('***********Pressao no fundo maior*******\n')
end
Qr = ((1).*((Entradas(3) - P_BHP).^(0.5))*P.k_resv(1));
% Vaz�o atrav�s da v�lvula choke.
%if(h0_anular > 0)
   %Qc = ((Entradas(2)).*((P_wh-Entradas(4)).^0.5)*P.k_choke(1));
%else
%    Qc = 0.5;
%end
% Equa��es diferenciais
dydt(1,1) = (Qr - q_media)/(P.Aanular*1000)/P.L_poco;
dydt(2,1) = ((P.B(2)/P.V(2))*(q_media-Qc)/10^3)./P.DeltaP_max;
dydt(3,1) = (1000*(P.med_Val(3)/(P.med_Val(4)*P.med_Val(1)))*(P_BHP*10^5 - P_wh*10^5 - sum(F) - ...
            (P.rho(1)*P.g*(P.h(1)-P.h(2)))-(P.rho(2)*P.g*(P.h(2)-P.h(3)))+ (dp)))/P.Delta_vazao;

  
end






        
        
        
