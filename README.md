# multiplex
Projet de communication avec méthode de multiporteuse  

## Réalisés

## RSB par sous-porteuse

### les courbes de TEB pour un canal AWGN

On fixe un RSB global du signal, qui va faire calculer la variance du bruit en fonction de celle du signal (( peut-être seulement implicitement : je suis pas sûr que ça aussi soit fait/ait besoin d'être))

Ca marche uniquement parce que le canal est constant en freq/égal à 1

Si canal Rayleigh, deux potentielles solutions (je pense)

Soit on s'assure de générer un canal de Rayleigh dont les coefficients sont aléatoires mais en s'assurant qu'il soit de moyenne nulle et avec variance normalisée
(Comme ce qu'on faisait en traitement d'image avec filtre 3D)

Soit on fixe le RSB du signal global, mais on connait pas le RSB effectifs des sous-porteuses passées dans le canal, qui va potentiellement réduire ce RSB des sous-porteuse. On va donc :

Fixer le RSB global pour la boucle for, mais afficher le RSB moyen des sous-porteuses, c'est-à-dire le RSB, mais après passage dans le canal et donc influence sur les sous-porteuses 
