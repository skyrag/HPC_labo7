#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>         

#define MAX_KMER 100

typedef struct KmerEntry {
    int count;
    struct KmerEntry *next;
    char kmer[];
} KmerEntry;

typedef struct {
    KmerEntry  **buckets;
    omp_lock_t  *locks;  // tableau de locks, un par bucket
    size_t       nbuckets;
    size_t       nentries;
    int          k;
} HashTable;

static long long now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long long)ts.tv_sec * 1000LL + (long long)ts.tv_nsec / 1000000LL;
}

static uint64_t hash_ptr(const char *s, int len) {
    uint64_t h = 5381;
    for (int i = 0; i < len; i++)
        h = ((h << 5) + h) + (unsigned char)s[i];
    return h;
}

static HashTable *hash_create(int k, size_t nbuckets) {
    HashTable *t = malloc(sizeof(HashTable));
    if (!t) exit(1);
    t->nbuckets = nbuckets;
    t->nentries = 0;
    t->k        = k;

    t->buckets = calloc(nbuckets, sizeof(KmerEntry *));
    if (!t->buckets) exit(1);

    t->locks = malloc(nbuckets * sizeof(omp_lock_t));
    if (!t->locks) exit(1);
    for (size_t i = 0; i < nbuckets; i++)
        omp_init_lock(&t->locks[i]);

    return t;
}

static void hash_free(HashTable *t) {
    if (!t) return;
    for (size_t i = 0; i < t->nbuckets; i++) {
        omp_destroy_lock(&t->locks[i]);   
        KmerEntry *e = t->buckets[i];
        while (e) {
            KmerEntry *n = e->next;
            free(e);
            e = n;
        }
    }
    free(t->locks);
    free(t->buckets);
    free(t);
}

static void hash_increment(HashTable *t, const char *ptr) {
    uint64_t h   = hash_ptr(ptr, t->k);
    size_t   idx = (size_t)(h & (t->nbuckets - 1));

    omp_set_lock(&t->locks[idx]);  // verrouillage du bucket

    KmerEntry *e = t->buckets[idx];
    while (e) {
        if (memcmp(e->kmer, ptr, t->k) == 0) {
            e->count++;
            omp_unset_lock(&t->locks[idx]);
            return;
        }
        e = e->next;
    }

    // in case not found
    KmerEntry *ne = malloc(sizeof(KmerEntry) + t->k + 1);
    if (!ne) { omp_unset_lock(&t->locks[idx]); exit(1); }
    memcpy(ne->kmer, ptr, t->k);
    ne->kmer[t->k] = '\0';
    ne->count  = 1;
    ne->next   = t->buckets[idx];
    t->buckets[idx] = ne;

    #pragma omp atomic
    t->nentries++;

    omp_unset_lock(&t->locks[idx]);  // déverrouillage du bucket
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <k>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int k = atoi(argv[2]);
    if (k <= 0) { fprintf(stderr, "k must be > 0\n"); return EXIT_FAILURE; }
    k = (k > MAX_KMER) ? MAX_KMER : k;

    FILE *file = fopen(argv[1], "rb");
    if (!file) { perror("fopen"); return EXIT_FAILURE; }
    if (fseek(file, 0, SEEK_END) != 0) { perror("fseek"); fclose(file); return EXIT_FAILURE; }
    long file_size = ftell(file);
    if (file_size < 0) { perror("ftell"); fclose(file); return EXIT_FAILURE; }
    if (file_size < k) { fclose(file); return 0; }
    rewind(file);

    char *buffer = malloc((size_t)file_size);
    if (!buffer) { perror("malloc"); fclose(file); return EXIT_FAILURE; }
    size_t nread = fread(buffer, 1, (size_t)file_size, file);
    if (nread != (size_t)file_size) file_size = (long)nread;
    fclose(file);

    size_t    nbuckets = 1u << 16;   
    HashTable *table   = hash_create(k, nbuckets);

    long nmers = file_size - k;      

    long long startTime = now_ms();

    #pragma omp parallel for schedule(dynamic, 512) default(none) \
            shared(buffer, table, nmers, k)
    for (long i = 0; i <= nmers; i++) {
        hash_increment(table, buffer + i);
    }

    long long endTime = now_ms();


    printf("Results:\n");
    for (size_t i = 0; i < table->nbuckets; i++) {
        for (KmerEntry *e = table->buckets[i]; e; e = e->next)
            printf("%s: %d\n", e->kmer, e->count);
    }

    printf("Time taken: %lld ms  |  threads: %d\n",
           endTime - startTime, omp_get_max_threads());

    hash_free(table);
    free(buffer);
    return 0;
}