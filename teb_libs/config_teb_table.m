%% Parametres de config du tableau

R = 1; % ici seuleument

EbN0dB_min  = 0; % Minimum de EbN0
EbN0dB_max  = 15; % Maximum de EbN0
EbN0dB_step = 1;% Pas de EbN0

nbr_erreur  = 200;  % Nombre d'erreurs à observer avant de calculer un BER
nbr_bit_max = 100e6;% Nombre de bits max à simuler
ber_min     = 3e-6; % BER min

EbN0dB = EbN0dB_min:EbN0dB_step:EbN0dB_max;     % Points de EbN0 en dB à simuler
EbN0   = 10.^(EbN0dB/10);% Points de EbN0 à simuler
EsN0   = R*log2(Mmod)*EbN0; % Points de EsN0
EsN0dB = 10*log10(EsN0); % Points de EsN0 en dB à simuler

% Là, on a configuré le EbN0dB_sousporteuse
EbN0dB_sousporteuse = EbN0dB;
EbN0dB = EbN0dB_sousporteuse./(abs(H).^2);


%% Variables init

debits_encod = [];
debits_decod = [];

%% Initialisation des vecteurs de r�sultats
ber = zeros(1,length(EbN0dB));
Pe = qfunc(sqrt(2*EbN0));

%% Préparation de l'affichage de la courbe

figure(1)

h_ber = semilogy(EbN0dB, ber,'XDataSource','EbN0dB', 'YDataSource','ber');
hold all
ylim([1e-6 1])
grid on
xlabel('$\frac{E_b}{N_0}$ en dB','Interpreter', 'latex', 'FontSize',14)
ylabel('TEB','Interpreter', 'latex', 'FontSize',14)

%% Préparation de l'affichage en console

msg_format = '|   %7.2f  |   %9d   |  %9d | %2.2e |  %8.2f kO/s |   %8.2f kO/s |   %8.2f s |\n';

fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')
msg_header =  '|  Eb/N0 dB  |    Bit nbr    |  Bit err   |   TEB    |    Debit Tx    |     Debit Rx    | Tps restant  |\n';
fprintf(msg_header);
fprintf(      '|------------|---------------|------------|----------|----------------|-----------------|--------------|\n')
