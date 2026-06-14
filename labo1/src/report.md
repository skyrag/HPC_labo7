
# **Rapport – Analyse automatique de signaux ECG**

## Introduction

L’objectif de ce laboratoire est de concevoir un système capable de détecter automatiquement les pics R dans des signaux ECG et de calculer les intervalles RR.

L’approche choisie repose sur une **version simplifiée de l’algorithme de Pan–Tompkins** : filtrage, dérivée, mise au carré, intégration sur fenêtre glissante, détection par seuil. Cette version est volontairement **non optimisée**, mais permet de valider le pipeline complet.

## Méthode

### Pipeline de traitement

1. **Suppression du DC** : retrait de la composante continue du signal.
2. **Filtrage passe-bande** : moyenne mobile puis filtre passe-haut pour isoler le complexe QRS.
3. **Dérivée du signal** : accentuation des pentes.
4. **Mise au carré** : amplification des variations rapides.
5. **Intégration sur fenêtre glissante** : moyenne sur 150 ms pour lisser le signal.
6. **Recentrement et gain** : remise à zéro du minimum et amplification du signal pour un seuil cohérent.
7. **Détection des pics R** : maximum local dépassant 50 % du pic global, avec période réfractaire de 200 ms.
8. **Calcul des intervalles RR** : différence en secondes entre pics consécutifs.

### Buffer et structure

* Deux buffers `current_sample1` et `current_sample2` sont utilisés pour stocker les étapes intermédiaires.
* Les pics R et les intervalles RR sont stockés dans des structures standard `ECG_Peaks` et `ECG_Intervals`.

## Performances

a minimum sur plusieurs itération le code a mis 0.000160 s a s'éxécuter

## Discussion

N'ayant jamais vraiment suivi de cours de C j'ai essayer de rattrapper une bonne partie de PRG2 pour faire ce labo mais cela m'a pris un certain temps et je n'ai pas eu le temps pour finir proprement et je m'en excuse, ce qu'y fait que mes points sont un peu décaler.


## Conclusion

Cette implémentation simple permet de :

* détecter les pics R sur des signaux ECG réels
* calculer les intervalles RR

