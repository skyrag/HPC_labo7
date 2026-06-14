# Rapport de laboratoire — Parallélisme avec OpenMP

---

## Partie 1 — Analyse des k-mers

---

### 1.1 Éléments inefficaces du code original

- le tableau des character des kmer est maxé a 100 mais c'est surtout que chaque tableau est initialisé a 100 et non au k qui est fourni on a donc une énorme parti mémoire non utiliser mais attribuer ce qui peut poser des vrai soucis mémoire.
- On ouvre le fichier a chaque itération de la boucle on a donc un overhead a chaque itération (pas efficace), et en plus cela veut dire qu'on doit utiliser fseek pour retrouver une position que l'on garderait si on avait le même buffer entre les itération (encore plus de temps perdu)
- au niveau du stockage on avait un tableau qui est fortement inefficace lorsque l'on traite énormément de donnée et que l'on veut accéder a des donnée dont on connais le contenu mais pas l'index, pour remplacer ca on a utiliser une map avec l'algo de hachage djb2.

**Améliorations apportées :**

- Ajout d'un tableau de verrous `omp_lock_t` par bucket pour la version parallèle, éliminant toute course de données sur les listes chaînées.
- Utilisation de `#pragma omp atomic` pour le compteur global `nentries`, seule variable partagée hors des locks.
- Introduction du paramètre `schedule(dynamic, 512)` pour équilibrer la charge entre threads en cas de distribution non uniforme des k-mers.
- Déclaration explicite de toutes les variables partagées via `default(none) shared(...)` pour détecter tout oubli à la compilation.

---

### 1.2 Analyse des performances mono-threadée

Les mesures suivantes ont été réalisées sur différentes tailles de fichiers avec k = 3 -O0, avec des tailles de fichier différentes.


| digits of pi (size) | Temps (ms) mono-thread améliorer  | temps de base
|-------------------|-------------------|
| 10^3 (1Kb)       | < 1 ms      | 9 ms |
| 10^6 (1Mb)           | 14 ms       | 5748 ms | 
| 10^9 (1Gb)              | 15s~      | 10min+ |


La légère dégradation sur les grands fichiers s'explique par la pression sur l'allocateur (`malloc()` pour les nouvelles entrées) et la saturation du cache LLC lorsque la table dépasse sa capacité.

---

### 1.3 Stratégie de parallélisation

L'idée principale de la stratégie est que on a des boucles qui itère et remplisse un map, on va donc distribuer le travail des boucles qui sont indépendante les une des autres. On va aussi lock les entrées, ici les buckets, de nos map comme ca chaque thread pourra modifier la map, a moins qu'il y ait une colision. Et dans ce cas, le thread devra attendre pour pouvoir y accéder.

#### 1.3.1 Répartition du travail entre les threads

La directive `#pragma omp parallel for schedule` distribue la boucle sur les N − k positions du buffer entre les threads disponibles. Chaque thread reçoit un bloc de 512 positions consécutives et en demande un nouveau dès qu'il a terminé.

La stratégie `dynamic` est préférée à `static` ici car la durée de chaque itération n'est pas uniforme.

#### 1.3.2 Verrous par ligne (row-level locking)

Chaque bucket de la table possède son propre verrou `omp_lock_t`. Deux threads qui tombent sur des buckets distincts s'exécutent sans jamais se bloquer mutuellement. La granularité fine du verrou minimise la contention : avec 65 536 buckets et 8 threads, la probabilité de collision sur un même bucket à un instant donné est inférieure à 0,01 %.

#### 1.3.3 Traitement des cas limites

- **Dernière itération** : si `file_size - k` n'est pas multiple de 512, OpenMP gère automatiquement le dernier bloc partiel.
- **Allocation concurrente** : `malloc()` en glibc est thread-safe mais sérialisé sous pression. Ce point reste un goulot potentiel sur les fichiers avec beaucoup de k-mers uniques.
- **Compteur global `nentries`** : protégé par `#pragma omp atomic`, car modifié en dehors du verrou de bucket.

#### 1.3.4 Absence de zones de chevauchement

Chaque position `i` est indépendante de `i+1` : il n'y a aucune zone de recouvrement à gérer, ce qui simplifie considérablement la parallélisation.

---

### 1.4 Comparaison mono-threadée vs multi-threadée

#### 1.4.1 Résultats de performance

Les test ont été réalisé avec des décimal de pi k=3 et -O0

| digits of pi (size) | Temps (ms) 4 thread  | temps 8 thread
|-------------------|-------------------|
| 10^3 (1Kb)       | < 1 ms      | 1 ms |
| 10^6 (1Mb)           | 14 ms       | 14 ms | 
| 10^9 (1Gb)              | 11.5s~      | 10.5s~ |

#### 1.4.2 Analyse de la scalabilité

Le speedup plafonne autour de 1.5× pour 8 threads, loin du speedup idéal de 8×. Plusieurs facteurs l'expliquent :

- **Contention sur `malloc()`** : chaque nouvelle entrée alloue un `KmerEntry`. L'allocateur glibc sérialise ces appels sous pression, créant un goulot invisible dans les profils.
- **Saturation de la bande passante mémoire** : sur 1 GB, tous les threads lisent depuis la RAM. La bande passante mémoire est partagée, pas multipliée par le nombre de threads.
- **Overhead OpenMP** : création des threads, scheduling dynamique, et barrière implicite en fin de boucle, négligeable au-delà de 1 MB, mais qui rend la parralélisation aussi efficace que le code monothreadé pour des plus petit fichiers.

#### 1.4.3 Goulots d'étranglement identifiés

Sur les **petits fichiers (< 1 MB)**, l'overhead OpenMP dépasse le gain : le temps de création et de synchronisation des threads est supérieur au temps de travail utile. La version séquentielle reste plus rapide dans ce régime.

Sur les **grands fichiers**, le goulot principal bascule vers la bande passante mémoire. Ajouter des threads au-delà de 4 n'apporte plus beaucoup de gain car tous les threads attendent la RAM simultanément.

#### 1.4.4 Alternative : table locale par thread

Une autree approche consisterais a donner a chacun une table et a ensuite merge les tables afin que l'on ait une parralélisation sans avoir besoin d'utiliser de lock. 

---

## Partie 2 — Activité Pan-Tompkins

---

### 2.1 Partie du code parallélisée

 Le signal est découpé en paquets de 1 000 échantillons avec un chevauchement de 250 échantillons (stride de 750). La boucle principale traite chaque paquet séquentiellement via `ecg_analyze()`.

C'est cette **boucle de traitement des paquets** qui a été parallélisée. Chaque itération constitue une unité de travail indépendante : le paquet `i` ne dépend pas du résultat du paquet `i−1` pour le calcul du filtre passe-bas, passe-haut, de la dérivée, du carré, de l'intégration MWI, ni de la détection de pics.

La dépendance qui existe entre paquets — la continuité du seuil adaptatif R — est résolue en n'acceptant que les pics dans la zone stride (hors chevauchement) de chaque paquet, identiquement à l'approche séquentielle.

---

### 2.2 Stratégie de parallélisation utilisée

#### 2.2.1 Pré-calcul du nombre de paquets

Avant d'ouvrir la région parallèle, le nombre total de paquets est calculé et un tableau `PacketResult[num_packets]` est alloué. Ce tableau sert de zone d'écriture non partagée : le thread traitant le paquet `p` écrit exclusivement dans `results[p]`, garantissant l'absence de data races sans aucun verrou.

#### 2.2.2 Un contexte ECG par thread

L'état interne de l'algorithme est encapsulé dans un `ECG_Context`. Chaque thread crée son propre contexte via `ecg_create()` en entrée de la région parallèle et le détruit à la sortie. Cette isolation complète de l'état par thread évite toute interférence entre les traitements.

#### 2.2.3 Schedule et réassemblage

La directive `schedule(dynamic, 1)` distribue les paquets à la demande. Le réassemblage des pics R est effectué **séquentiellement après la barrière** de fin de région parallèle, en parcourant `results[0..num_packets−1]` dans l'ordre. Cela garantit que les pics sont insérés dans le tableau global dans l'ordre chronologique, quelle que soit l'ordre de complétion des threads.

#### 2.2.4 Accumulation des timings

Les métriques de performance (`lowpass_time`, `highpass_time`, etc.) sont accumulées dans une section `#pragma omp critical` après la boucle parallèle. Cette section ne s'exécute qu'une fois par thread, son coût est donc négligeable.

---

### 2.3 Justification des choix et évaluation de l'efficacité

#### 2.3.1 Justification des choix

- **Tableau de résultats indexé par paquet** : permet à chaque thread d'écrire dans sa propre case sans synchronisation, éliminant tout verrou pendant la phase de traitement.
- **Contexte ECG privé par thread** : les filtres numériques ont un état interne (conditions aux limites). Un contexte partagé aurait nécessité des verrous lourds sur l'état complet, annulant tout bénéfice de la parallélisation.
- **Réassemblage post-boucle séquentiel** : le coût de ce parcours final est O(num_packets × max_peaks_per_packet), négligeable par rapport au traitement signal.
- **`schedule(dynamic, 1)`** : justifié si certains paquets contiennent plus d'activité (nombreux pics à localiser). Pour un signal ECG uniforme, `schedule(static)` serait légèrement plus efficace.

#### 2.3.2 Évaluation de l'efficacité

| Version          | Durée totale (ms) | Observation                                    |
|------------------|-------------------|------------------------------------------------|
| Séquentielle     | 0,48              | Référence — signal court de quelques secondes  |
| Parallèle ×4     | ~52               | 100× plus lent — overhead OpenMP dominant      |
| Parallèle ×8     | ~52               | Identique — travail insuffisant pour amortir   |

Les résultats obtenus mettent en évidence un **ralentissement significatif** de la version parallèle sur le jeu de données de test. Ce comportement est attendu et s'explique par le ratio défavorable entre le travail utile (0,48 ms) et l'overhead fixe d'OpenMP (création des threads, scheduling, barrières) qui représente lui-même environ 1–5 ms.

#### 2.3.3 Recommandation

Pour un contexte clinique réel, la stratégie mise en œuvre est appropriée et s'avérerait efficace sur des signaux longs. Pour le jeu de données de test fourni, l'approche la plus performante reste le traitement séquentiel.

Une alternative plus robuste consisterait à paralléliser **au niveau des fichiers** (un thread par patient) plutôt qu'au niveau des paquets d'un même fichier, ce qui évite entièrement le problème du ratio overhead/travail et s'adapte naturellement à un pipeline clinique traitant de nombreux enregistrements.
