% Método de Hardy-Cross para resolução de redes malhadas
clear all;
clc;

% Entrada
massaespec = 1000;               % Massa específica (kg/m^3)
viscosidade = 1.002 * 10^(-3);   % Viscosidade (Pa.s)
iteracoes = 0;

% Comprimento (m), diâmetro (m), rugosidade (mm), vazão (m^3/s)
aneis(:,:,1) = [
    1000, 1, 0.0015, 1;  % 1
    925,  0.75,  0.0015, 6;   % 2
    1000,  0.75,  0.0015, 3;  % 3
    925,  0.75,  0.0015, -1; % 4
    1,    1,    1,     1;
    1,    1,    1,     1        % Extra linha para completar a matriz
];

aneis(:,:,2) = [
   925, 0.75, 0.0015, -6;  % 2
   1000, 0.75, 0.0015, -1;  % 9
   1000, 1, 0.0015, -6;  % 10
   925, 1, 0.0015, 6;  % 11
    1,    1,    1,     1; % Extra linha para completar a matriz
    1,    1,    1,     1
];

aneis(:,:,3) = [
  1000, 0.75, 0.0015,-3 ;  % 3
   350, 0.5, 0.0015,-3 ;  % 5
   671, 0.5, 0.0015, -2;  % 6
   400, 0.5, 0.0015, -1;  % 7
   650, 0.5, 0.0015, 1;  % 8
   1,    1,    1,     1 % Extra linha para completar a matriz
];

aneis(:,:,4) = [
   650, 0.5, 0.0015,-1 ;  % 8
   1000, 0.75, 0.0015,1 ;  % 9
   800, 1, 0.0015,6 ;  % 12
   650, 1, 0.0015, 4;  % 14
   800, 0.5, 0.0015, 2;  % 21
   1000, 0.5, 0.0015, -1;  % 22
];

aneis(:,:,5) = [
   763,0.75 , 0.0015, 1;  % 13
   650, 1, 0.0015, -4;  % 14
   400, 0.75, 0.0015, -0.5;  % 15
    1,    1,    1,     1; % Extra linha para completar a matriz
    1,    1,    1,     1;   
    1,    1,    1,     1
    ];

aneis(:,:,6) = [
   125, 0.5, 0.0015,0.5 ;  % 18
   800, 0.5, 0.0015,1 ;  % 19
   125, 0.5, 0.0015,-2 ;  % 20
   800, 0.5, 0.0015, -2;  % 21
   1,    1,    1,     1; % Extra linha para completar a matriz
   1,    1,    1,     1     
];

aneis(:,:,7) = [
   400, 0.75, 0.0015, 0.5;  % 15
   125, 0.75, 0.0015, 0.5;  % 16
   400, 0.75, 0.0015, 0.5;  % 17
   125, 0.5, 0.0015, -0.5;  % 18
   1,    1,    1,     1; % Extra linha para completar a matriz
   1,    1,    1,     1       
];



% Número de trechos dos anéis
numT = [4, 4, 5, 6, 3, 4,4];
% Número de anéis
numA = 7;

% Inicializando variáveis
dQ = zeros(1, numA);
H = zeros(numA, max(numT));

% Para armazenar os valores de delta Q em cada iteração
deltaQ1_vals = [];
deltaQ2_vals = [];
deltaQ3_vals = [];
deltaQ4_vals = [];
deltaQ5_vals = [];
deltaQ6_vals = [];
deltaQ7_vals = [];

while 1
    % Trechos repetidos: (6)
    aneis(2,4,1) = aneis(2,4,1) - dQ(2);%2
    aneis(1,4,2) = aneis(5,4,2) - dQ(1);
    aneis(3,4,1) = aneis(3,4,1) - dQ(3);%3
    aneis(1,4,3) = aneis(1,4,3) - dQ(1);
    aneis(5,4,3) = aneis(5,4,3) - dQ(4);%8
    aneis(1,4,4) = aneis(1,4,4) - dQ(3);
    aneis(2,4,2) = aneis(2,4,2) - dQ(4);%9
    aneis(2,4,4) = aneis(2,4,4) - dQ(2);
    aneis(5,4,4) = aneis(5,4,4) - dQ(6);%21
    aneis(4,4,6) = aneis(4,4,6) - dQ(4);
    aneis(4,4,4) = aneis(4,4,4) - dQ(5);%14
    aneis(2,4,4) = aneis(2,4,4) - dQ(4);
    aneis(1,4,6) = aneis(1,4,6) - dQ(7);%18
    aneis(4,4,7) = aneis(4,4,7) - dQ(6);
    aneis(3,4,5) = aneis(3,4,5) - dQ(7);%15
    aneis(1,4,7) = aneis(1,4,7) - dQ(5);


    for j = 1:numA
        aneis(:,4,j) = aneis(:,4,j) + dQ(j);
        
        for i = 1:numT(j)
            H(j,i) = perdacarga(massaespec, viscosidade, aneis(i,1,j), aneis(i,2,j), aneis(i,3,j), aneis(i,4,j));
        end
        
        HsobreQ(j,:) = H(j,:) ./ aneis(:,4,j)';
        dQ(j) = -sum(H(j,:)) / (2 * sum(HsobreQ(j,:)));
        somaH(j) = sum(H(j,:));
    end
    
    iteracoes = iteracoes + 1;
    
    % Armazena valores de delta Q em mililitros (m³/s * 10^-3)
    deltaQ1_vals(iteracoes) = dQ(1) * 10^3;
    deltaQ2_vals(iteracoes) = dQ(2) * 10^3;
    deltaQ3_vals(iteracoes) = dQ(3) * 10^3;
    deltaQ4_vals(iteracoes) = dQ(4) * 10^3;
    deltaQ5_vals(iteracoes) = dQ(5) * 10^3;
    deltaQ6_vals(iteracoes) = dQ(6) * 10^3;
    deltaQ7_vals(iteracoes) = dQ(7) * 10^3;
    
    % Condição de parada
    if abs(max(somaH)) < 1 * 10^(-5)
        break
    end
end

% Exibe os valores armazenados de ∆Q1 a ∆Q7 com 9 casas decimais
fprintf('Iteração   ∆Q1 (m³/s)×10-³       ∆Q2 (m³/s)×10-³       ∆Q3 (m³/s)×10-³       ∆Q4 (m³/s)×10-³       ∆Q5 (m³/s)×10-³       ∆Q6 (m³/s)×10-³       ∆Q7 (m³/s)×10-³\n');
fprintf('-------------------------------------------------------------------------------------------------------------------------\n');
for iter = 1:iteracoes
    fprintf('%4d      %+15.9f     %+15.9f     %+15.9f     %+15.9f     %+15.9f     %+15.9f     %+15.9f\n', ...
        iter, deltaQ1_vals(iter), deltaQ2_vals(iter), deltaQ3_vals(iter), deltaQ4_vals(iter), deltaQ5_vals(iter), deltaQ6_vals(iter), deltaQ7_vals(iter));
end


fprintf('As vazões obtidas foram (desconsidere valores das linhas extras):\n');
aneis(:,4,:)

%Perda de Carga
function H = perdacarga(massaespec, viscosidade, comprimento, diametro, rugosidade, vazao)

    % Área da seção transversal do tubo
    area = pi * diametro^2 / 4;
    
    % Cálculo da velocidade do fluido
    velocidade = vazao / area;
    
    % Cálculo do número de Reynolds (Re)
    Re = massaespec * abs(velocidade) * diametro / viscosidade;
    
    % Equação de Swamee-Jain para o cálculo do fator de atrito
    fator = 0.25 / (log10(rugosidade / (3.7 * diametro * 10^3) + 5.74 / Re^0.9))^2;
    
    % Equação de Darcy-Weisbach para a perda de carga
    H = (velocidade / abs(velocidade)) * comprimento * fator * velocidade^2 / (9.81 * 2 * diametro);

end

