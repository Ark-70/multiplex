function bruit = calculerBruit(RSB, signal)
    % Calcul le bruit en fonction du RSB

    bruit_initial = randn(length(echantillon),1);


    puissance_signal = sum(real(echantillon).^2)/length(real(echantillon));
    puissance_bruit = sum(bruit_initial.^2)/length(bruit_initial);
    alpha = sqrt((puissance_signal/puissance_bruit) * 10^(-RSB/10));

    bruit = alpha * bruit_initial;
end
