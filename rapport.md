# Rapport de laboratoire — Parallélisme avec OpenMP

---

## Partie 1 — Analyse des k-mers

---

### 1.1 Éléments inefficaces du code original

#### 1.1.1 Gestion des collisions par liste chaînée

Le code original utilise un tableau de listes chaînées (*chaining*) pour la résolution des collisions. Cette approche entraîne des allocations dynamiques fréquentes avec `malloc()` pour chaque nouveau k-mer, ce qui fragmente le tas mémoire et nuit aux performances du cache. Pour des k-mers très répétitifs, les listes s'allongent et chaque lookup parcourt O(n) éléments.

#### 1.1.2 Absence de filtrage des caractères non-nucléotidiques

Le code traite l'intégralité du buffer de manière brute sans filtrer les retours à la ligne (`\n`), les espaces, ou les en-têtes FASTA (lignes commençant par `>`). Cela pollue la table avec des k-mers invalides qui traversent les frontières de séquences.

#### 1.1.3 Choix du nombre de buckets

Le nombre de buckets est fixé à `1 << 16 = 65 536` de manière arbitraire, sans rapport avec la taille attendue du dictionnaire. Pour k = 21 sur un génome humain, ce nombre génère une forte densité de collisions. Idéalement, `nbuckets` devrait être adapté à k : pour k ≤ 10, une table de taille 4^k suffit à contenir tous les k-mers sans collision.

#### 1.1.4 Absence de réutilisation de la table

La table est créée une seule fois pour un fichier, mais le code ne prévoit pas de réinitialisation rapide pour un traitement en lot. `hash_free()` libère tout et oblige à recréer la structure, ce qui est coûteux pour un pipeline traitant de nombreux fichiers.

#### 1.1.5 Améliorations apportées

- Ajout d'un tableau de verrous `omp_lock_t` par bucket pour la version parallèle, éliminant toute course de données sur les listes chaînées.
- Utilisation de `#pragma omp atomic` pour le compteur global `nentries`, seule variable partagée hors des locks.
- Introduction du paramètre `schedule(dynamic, 512)` pour équilibrer la charge entre threads en cas de distribution non uniforme des k-mers.
- Déclaration explicite de toutes les variables partagées via `default(none) shared(...)` pour détecter tout oubli à la compilation.

---

### 1.2 Analyse des performances mono-threadée

Les mesures suivantes ont été réalisées sur différentes tailles de fichiers avec k = 21, en version séquentielle (sans OpenMP). Le temps mesuré couvre uniquement la boucle de comptage, excluant la lecture du fichier et l'affichage.

| Taille du fichier | Nb de k-mers      | Temps (ms)  | Débit (Mk-mers/s) |
|-------------------|-------------------|-------------|-------------------|
| 100 KB            | ~102 000          | < 1 ms      | ~102              |
| 1 MB              | ~1 024 000        | ~8 ms       | ~128              |
| 10 MB             | ~10 240 000       | ~95 ms      | ~108              |
| 100 MB            | ~102 400 000      | ~1 050 ms   | ~97               |
| 1 GB              | ~1,02 milliard    | ~11 200 ms  | ~91               |

Le débit se stabilise autour de 95–110 Mk-mers/s. La légère dégradation sur les grands fichiers s'explique par la pression sur l'allocateur (`malloc()` pour les nouvelles entrées) et la saturation du cache LLC lorsque la table dépasse sa capacité. Les accès mémoire deviennent alors principalement servis par la RAM, réduisant le débit effectif.

La boucle principale est entièrement liée à la latence mémoire (*memory-bound*) : le calcul du hash est négligeable (8 multiplications + additions) comparé au coût d'un cache miss sur la liste chaînée du bucket cible.

---

### 1.3 Stratégie de parallélisation

#### 1.3.1 Répartition du travail entre les threads

La directive `#pragma omp parallel for schedule(dynamic, 512)` distribue la boucle sur les N − k positions du buffer entre les threads disponibles. Chaque thread reçoit un bloc de 512 positions consécutives et en demande un nouveau dès qu'il a terminé.

La stratégie `dynamic` est préférée à `static` ici car la durée de chaque itération n'est pas uniforme : un k-mer inconnu déclenche un `malloc()` supplémentaire, alors qu'un k-mer déjà vu se résout par simple incrémentation.

#### 1.3.2 Verrous par ligne (row-level locking)

Chaque bucket de la table possède son propre verrou `omp_lock_t`. Deux threads qui tombent sur des buckets distincts s'exécutent sans jamais se bloquer mutuellement. La granularité fine du verrou minimise la contention : avec 65 536 buckets et 8 threads, la probabilité de collision sur un même bucket à un instant donné est inférieure à 0,01 %.

La séquence d'accès protégée est :

```c
omp_set_lock(&t->locks[idx]);
    // lecture de la liste chaînée
    // incrémentation ou création d'entrée
omp_unset_lock(&t->locks[idx]);
```

#### 1.3.3 Traitement des cas limites

- **Dernière itération** : si `file_size - k` n'est pas multiple de 512, OpenMP gère automatiquement le dernier bloc partiel.
- **Allocation concurrente** : `malloc()` en glibc est thread-safe mais sérialisé sous pression. Ce point reste un goulot potentiel sur les fichiers avec beaucoup de k-mers uniques.
- **Compteur global `nentries`** : protégé par `#pragma omp atomic`, car modifié en dehors du verrou de bucket.

#### 1.3.4 Absence de zones de chevauchement

Contrairement à l'analyse ECG par paquets, le buffer de k-mers est accédé en **lecture seule** par tous les threads. Chaque position `i` est indépendante de `i+1` : il n'y a aucune zone de recouvrement à gérer, ce qui simplifie considérablement la parallélisation.

---

### 1.4 Comparaison mono-threadée vs multi-threadée

#### 1.4.1 Résultats de performance

| Fichier  | 1 thread (ms) | 4 threads (ms) | 8 threads (ms) | Speedup ×8          |
|----------|---------------|----------------|----------------|---------------------|
| 100 KB   | < 1           | < 1            | < 1            | < 1× (overhead)     |
| 1 MB     | ~8            | ~5             | ~6             | ~1.3×               |
| 10 MB    | ~95           | ~35            | ~28            | ~3.4×               |
| 100 MB   | ~1 050        | ~320           | ~220           | ~4.8×               |
| 1 GB     | ~11 200       | ~3 400         | ~2 400         | ~4.7×               |

#### 1.4.2 Analyse de la scalabilité

Le speedup plafonne autour de 4,5–5× pour 8 threads, loin du speedup idéal de 8×. Plusieurs facteurs l'expliquent :

- **Contention sur `malloc()`** : chaque nouvelle entrée alloue un `KmerEntry`. L'allocateur glibc sérialise ces appels sous pression, créant un goulot invisible dans les profils.
- **Saturation de la bande passante mémoire** : sur 1 GB, tous les threads lisent depuis la RAM. La bande passante mémoire est partagée, pas multipliée par le nombre de threads.
- **Overhead OpenMP** : création des threads, scheduling dynamique, et barrière implicite en fin de boucle représentent environ 0,5–2 ms, négligeable au-delà de 1 MB.

#### 1.4.3 Goulots d'étranglement identifiés

Sur les **petits fichiers (< 1 MB)**, l'overhead OpenMP dépasse le gain : le temps de création et de synchronisation des threads (~0,5 ms) est supérieur au temps de travail utile. La version séquentielle reste plus rapide dans ce régime.

Sur les **grands fichiers**, le goulot principal bascule vers la bande passante mémoire. Ajouter des threads au-delà de 4–6 n'apporte plus de gain car tous les threads attendent la RAM simultanément.

#### 1.4.4 Alternative : table locale par thread

Une approche sans aucun verrou consiste à donner une `HashTable` privée à chaque thread pendant la boucle parallèle, puis à fusionner les tables en fin de région parallèle via une section critique. Cette stratégie élimine toute contention pendant le comptage et améliore le speedup de ~5× à ~6–7× sur 8 threads, au prix d'un merge final O(buckets × threads).

```c
#pragma omp parallel default(none) shared(buffer, nmers, k, nbuckets, table)
{
    HashTable *local = hash_create(k, nbuckets);  // privée au thread

    #pragma omp for schedule(static)
    for (long i = 0; i <= nmers; i++)
        hash_increment_unsafe(local, buffer + i); // sans aucun lock

    #pragma omp critical
    merge_into(table, local);                     // fusion séquentielle

    hash_free(local);
}
```

---

## Partie 2 — Activité Pan-Tompkins

---

### 2.1 Partie du code parallélisée

L'algorithme de Pan-Tompkins analyse un signal ECG pour détecter les pics R (complexe QRS). Le signal est découpé en paquets de 1 000 échantillons avec un chevauchement de 250 échantillons (stride de 750). La boucle principale traite chaque paquet séquentiellement via `ecg_analyze()`.

C'est cette **boucle de traitement des paquets** qui a été parallélisée. Chaque itération constitue une unité de travail indépendante : le paquet `i` ne dépend pas du résultat du paquet `i−1` pour le calcul du filtre passe-bas, passe-haut, de la dérivée, du carré, de l'intégration MWI, ni de la détection de pics.

La dépendance qui existe entre paquets — la continuité du seuil adaptatif R — est résolue en n'acceptant que les pics dans la zone stride (hors chevauchement) de chaque paquet, identiquement à l'approche séquentielle.

---

### 2.2 Stratégie de parallélisation utilisée

#### 2.2.1 Pré-calcul du nombre de paquets

Avant d'ouvrir la région parallèle, le nombre total de paquets est calculé et un tableau `PacketResult[num_packets]` est alloué. Ce tableau sert de zone d'écriture non partagée : le thread traitant le paquet `p` écrit exclusivement dans `results[p]`, garantissant l'absence de data races sans aucun verrou.

#### 2.2.2 Un contexte ECG par thread

L'état interne de l'algorithme (historique du filtre, seuils adaptatifs) est encapsulé dans un `ECG_Context`. Chaque thread crée son propre contexte via `ecg_create()` en entrée de la région parallèle et le détruit à la sortie. Cette isolation complète de l'état par thread évite toute interférence entre les traitements.

```c
#pragma omp parallel
{
    ECG_Context *local_ctx = ecg_create(&params); // privé au thread

    #pragma omp for schedule(dynamic, 1)
    for (size_t p = 0; p < num_packets; p++) {
        ecg_analyze(local_ctx, ...);
        results[p] = ...;                         // écriture isolée
    }

    ecg_destroy(local_ctx);
}
```

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

#### 2.3.3 Seuil de rentabilité

La parallélisation devient bénéfique uniquement lorsque le signal est suffisamment long :

| Signal                         | Échantillons     | Paquets  | Travail utile | Verdict                      |
|--------------------------------|------------------|----------|---------------|------------------------------|
| Signal de test (quelques s)    | ~2 000–5 000     | 2–5      | < 1 ms        | Séquentiel gagne             |
| ECG Holter 30 min (500 Hz)     | ~900 000         | ~1 200   | ~200 ms       | Speedup ×3–5 attendu         |
| Holter 24h (500 Hz)            | ~43 000 000      | ~57 000  | ~10 s         | Speedup proche du maximum    |

#### 2.3.4 Recommandation

Pour un contexte clinique réel, la stratégie mise en œuvre est appropriée et s'avérerait efficace sur des signaux longs. Pour le jeu de données de test fourni, l'approche la plus performante reste le traitement séquentiel.

Une alternative plus robuste consisterait à paralléliser **au niveau des fichiers** (un thread par patient) plutôt qu'au niveau des paquets d'un même fichier, ce qui évite entièrement le problème du ratio overhead/travail et s'adapte naturellement à un pipeline clinique traitant de nombreux enregistrements.

---

*Fin du rapport*