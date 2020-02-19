clear, close all, clc;

%% PARAMÊTRES
% -------------------------------------------------------------------------

% Paramêtres de la chaîne de com

PREFIX_CYCL_ON = 1;
BRUIT_ON = 0;
CANAL_TYPE = 'Rayleigh'; % 'Rayleigh' ou 'AWGN'
ANNULATION_ON = 0;
EGALISEUR_ON = 1;

% Constantes générales

M=2; % M-PSK
Ts = 0.05e-6;
Fe = 1/Ts;
RSB = 5; % Définit l'amplitude du bruit

% Constantes OFDM

K = 500; % symboles OFDM d'une trame OFDM
N = 128; % Nombre de sous-porteuses totales
garde = 16; % intervalle de garde
annulation = 4;
L = 50; % Composantes cheloues du filtre (si 2 => 2 dirac -> un cos)

%% PARAMÊTRES CALCULÉS
% -------------------------------------------------------------------------
% NbTrames = 2;
Nbits = N*K;
nbMod = Nbits/N;
if(garde < L)
    fprintf("Warning : De l'IES va apparaitre à cause d'un intervalle de garde trop court.\n");
end

% Blocs de com

PSKMod   = comm.PSKModulator(M,'BitInput',true,'PhaseOffset',0);
PSKDemod = comm.PSKDemodulator(M);
errorRate = comm.ErrorRate;

%% CHAÎNE DE COM
% -------------------------------------------------------------------------

%% ÉMETTEUR
% -------------------------------------------------------------------------


bits = randi([0 1],Nbits,1);

%% modulation

symboles = step(PSKMod, bits);

%% IFFT
% (à faire sur matrice)

echantillon = zeros(length(symboles),1);
for i=1:nbMod
    echantillon(N*(i-1)+1:N*i) = ifft(symboles(N*(i-1)+1:N*i))*sqrt(N);

    % Intervalle d'annulation vers fe/2 et -fe/2 (les trucs sans intérêt)
    if (ANNULATION_ON)
        for g=-annulation+1:annulation
            echantillon(N*(i-1)+1+(N/2)+g) = 0;
        end
    end
end


%% Préfix
matriceTrames = reshape(echantillon, [N K]);
test1 = matriceTrames;
if (PREFIX_CYCL_ON)

    % Série -> Parallèle
    matriceTrames = reshape(echantillon, [N K]); % premiere trame = 1ere colonne

    suffixACopier = matriceTrames(end-garde+1:end,:);
    matriceTrames = [ suffixACopier ; matriceTrames ];

    % Parallèle -> Série
    echantillon = reshape(matriceTrames, (garde+N)*K, 1);
else
    garde = 0;
end


%% CANAL
% -------------------------------------------------------------------------

if (CANAL_TYPE == 'Rayleigh')
    % Génération de gaussiennes complexes comme composantes de canal
    h = sqrt(1/2*L)*(randn(1,L)+1j*randn(1,L));
elseif (CANAL_TYPE == 'AWGN')
    h = 1;
end

H = fft(h, N);

if (BRUIT_ON)
    bruit = calculerBruit(RSB, echantillon);
    y = echantillon + bruit;
else
    y = echantillon;
end

y = filter(h, 1, y);

%% RÉCEPTEUR
% -------------------------------------------------------------------------

% Série -> Parallèle

matriceTrames = reshape(y, [N+garde K]); % premiere trame = 1ere colonne

%% On enlève le prefix

matriceTrames = matriceTrames(garde+1:end,:);
test2 = matriceTrames;
% y = reshape(matriceTrames, N*K, 1);

%% FFT

% On estime qu'on connait h

% je prend la tf de h sur (n) points, pour trouver les différentes
% multiplications sur les sous-porteuses

%(censé être fait pas avec une boucle for mais sur une matrice)

matriceTrames = fft(matriceTrames, N)/sqrt(N);

%% Parallèle -> Série
symbolesRecus = reshape(matriceTrames, N*K, 1);

%% Egaliseur

%
if(EGALISEUR_ON)
    Hprep = repmat(H.', 500, 1); % .' parce qu'on veut pas le conjugué mais la transposée
    symbolesRecus = symbolesRecus./Hprep;
end

%% Repassage en bitstream

% Demodulation
bitstream = step(PSKDemod, symbolesRecus);
errorStats = errorRate(bits,bitstream);
fprintf('Error rate = %f\nNumber of errors = %d\n', ...
    errorStats(1), errorStats(2));


%% questions :

Tsp = N*Ts; %en s
df = 1/Tsp; %en Hz
%si Nu = N, la bande vaut df*N
bande = df*N;
Ds = Fe;

%% affichage
figure;
plot(real(echantillon));
figure;
subplot(2,1,1);
plot(real(echantillon(N*(1-1)+1:N*1)));
subplot(2,1,2);
plot(imag(echantillon(N*(1-1)+1:N*1)));

%% questions 2
figure;
subplot(2,1,1);
histogram(real(echantillon),100);
subplot(2,1,2);
histogram(imag(echantillon),100);


scatterplot(symbolesRecus);

rCorr = xcorr(real(echantillon));
iCorr = xcorr(imag(echantillon));
xCorr = xcorr(real(echantillon),imag(echantillon),'unbiased');
figure;
subplot(3,1,1);
plot(rCorr);
subplot(3,1,2);
plot(iCorr);
subplot(3,1,3);
plot(xCorr);

figure;
pwelch(echantillon,ones(1,N),0,4*N);
