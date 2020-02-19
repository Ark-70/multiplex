clear;
close all;
clc;
addpath("teb_libs");
%% PARAMÊTRES
% -------------------------------------------------------------------------

% Paramêtres de la chaîne de com

PREFIX_CYCL_ON = 1;
BRUIT_ON = 0;
CANAL_TYPE = 'Rayleigh'; % 'Rayleigh' ou 'AWGN'
ANNULATION_ON = 0;
EGALISEUR_ON = 1;
TEB_LOOP_ON = 1;
RSB_SOUS_PORTEUSE_EN_PLOT = 1;

% Constantes générales

Mmod=2; % M-PSK
Ts = 0.05e-6;
Fe = 1/Ts;
RSB = 5; % Définit l'amplitude du bruit

% Constantes OFDM

K = 500; % symboles OFDM d'une trame OFDM
N = 128; % Nombre de sous-porteuses totales
garde = 8; % intervalle de garde
annulation = 4;
L = 50; % Composantes cheloues du filtre (si 2 => 2 dirac -> un cos)
nbTrames = 64; % = Nbits/500/128 == toute la matrice temps-frequence OFDM

%% PARAMÊTRES CALCULÉS
% -------------------------------------------------------------------------

Nbits = N*K*nbTrames;
nbMod = Nbits/N;

if (garde < L)
    fprintf("Warning : De l'IES va apparaitre à cause d'un intervalle de garde trop court.\n");
end

if (~PREFIX_CYCL_ON)
    garde = 0;
end

if (~TEB_LOOP_ON)
    EbN0dB = RSB
end

% Blocs de com

PSKMod   = comm.PSKModulator(Mmod,'BitInput',true,'PhaseOffset',0);
PSKDemod = comm.PSKDemodulator(Mmod);
stat_erreur = comm.ErrorRate;

if (strcmp(CANAL_TYPE, 'Rayleigh')==0 )
    % Génération de gaussiennes complexes comme composantes de canal
    h = sqrt(1/2*L)*(randn(1,L)+1j*randn(1,L));
elseif (strcmp(CANAL_TYPE, 'AWGN')==0 )
    h = 1;
end
H = fft(h, N);

%% LOOP DE TEB/Variation EdBN0

config_teb_table;

for i_snr = 1:length(EbN0dB)
    RSB = EbN0dB(i_snr);
    reverseStr = ''; % Pour affichage en console

    stat_erreur.reset; % reset du compteur d'erreur
    err_stat    = [0 0 0]; % vecteur r�sultat de stat_erreur

    % demod_psk.Variance = awgn_channel.Variance;

    n_frame = 0;
    T_rx = 0;
    T_tx = 0;
    general_tic = tic;
    while (err_stat(2) < nbr_erreur && err_stat(3) < nbr_bit_max)
        n_frame = n_frame + 1;

%% CHAÎNE DE COM
% -------------------------------------------------------------------------
%% ÉMETTEUR
% -------------------------------------------------------------------------
tx_tic    = tic; % Mesure du débit d'encodage au début de l'emetteur

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

T_tx      = T_tx+toc(tx_tic); % Mesure du débit d'encodage fin de l'émetteur


%% CANAL
% -------------------------------------------------------------------------

if (BRUIT_ON)
    bruit = calculerBruit(RSB, echantillon);
    y = echantillon + bruit;
else
    y = echantillon;
end

y = filter(h, 1, y);

%% RÉCEPTEUR
% -------------------------------------------------------------------------
rx_tic  = tic; % Temps t début récepteur

% Série -> Parallèle

matriceTrames = reshape(y, [N+garde K*nbTrames]); % premiere trame = 1ere colonne

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
symbolesRecus = reshape(matriceTrames, N*K*nbTrames, 1);

%% Egaliseur

%
if(EGALISEUR_ON)
    Hprep = repmat(H.', K*nbTrames, 1); % .' parce qu'on veut pas le conjugué mais la transposée
    symbolesRecus = symbolesRecus./Hprep;
end

%% Repassage en bitstream

% Demodulation
bitstream = step(PSKDemod, symbolesRecus);
T_rx    = T_rx + toc(rx_tic);  % temps t à la fin du recepteur (pour le calcul du débit)
errorStats = stat_erreur(bits,bitstream);

% if(~TEB_LOOP_ON)
%     fprintf('Error rate = %f\nNumber of errors = %d\n', ...
%     errorStats(1), errorStats(2));
% end

%% Affichage d'une ligne de TEB
    err_stat   = step(stat_erreur, bits, bitstream); % Comptage des erreurs binaires

    if mod(n_frame,100) == 1

        debit_encod = err_stat(3)/8/T_tx/1e3; % 8 = octet?
        debit_decod = err_stat(3)/8/T_rx/1e3;

        display_str = sprintf(msg_format,...
            EbN0dB(i_snr),         ... % EbN0 en dB
            err_stat(3),           ... % Nombre de bits envoy�s
            err_stat(2),           ... % Nombre d'erreurs observ�es
            err_stat(1),           ... % BER
            debit_encod,... % D�bit d'encodage
            debit_decod,... % D�bit de d�codage
            toc(general_tic)*(nbr_erreur - min(err_stat(2),nbr_erreur))/nbr_erreur); % Temps restant
        fprintf(reverseStr);
        msg_sz =  fprintf(display_str);
        reverseStr = repmat(sprintf('\b'), 1, msg_sz);
    end
%% Fin de la Loop_de_TEB

        if(~TEB_LOOP_ON)
            break;
        end
    end%while

    debit_encod = err_stat(3)/8/T_tx/1e3;
    debit_decod = err_stat(3)/8/T_rx/1e3;

    debits_encod = [debits_encod debit_encod];
    debits_decod = [debits_decod debit_decod];

    display_str = sprintf(msg_format,EbN0dB(i_snr), err_stat(3), err_stat(2), err_stat(1), err_stat(3)/8/T_tx/1e3, err_stat(3)/8/T_rx/1e3, toc(general_tic)*(100 - min(err_stat(2),100))/100);
    fprintf(reverseStr);
    msg_sz =  fprintf(display_str);
    reverseStr = repmat(sprintf('\b'), 1, msg_sz);

    ber(i_snr) = err_stat(1);
    refreshdata(h_ber);
    drawnow limitrate

    if err_stat(1) < ber_min
        break
    end

end%for
%% Affichage final
load('BPSK_TEB_0_1_15.mat');
figure(1)
hold all
% même si pour le RSB global, c'est toujours un canal avec tel RSB global
% Mais pour ce graph on doit faire évoluer le RSB par sous-porteuse
% EbN0dB_sousporteuse = EbN0dB_sousporteuse*(abs(H)^2);
semilogy(EbN0dB,BPSKRayleighL16.data{1, 2});

% semilogy(EbN0dB,ber);
xlim([0 15])
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)
title('TEB pour chaque code');
legend(['OFDM sur canal Rayleigh L=16' newline '(EbN0 moyen des sous-porteuses)'], 'BPSK sur canal Rayleigh L=16');


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
