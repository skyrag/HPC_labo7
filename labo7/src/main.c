#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>




#define MAX_KMER 100

typedef struct KmerEntry {
    int count;
    struct KmerEntry *next;
    char kmer[];
} KmerEntry;

typedef struct {
    KmerEntry **buckets;
    size_t nbuckets;
    size_t nentries;
    int k;
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
    t->buckets = calloc(nbuckets, sizeof(KmerEntry*));
    if (!t->buckets) exit(1);
    t->nentries = 0;
    t->k = k;
    return t;
}

static void hash_free(HashTable *t) {
    if (!t) return;
    for (size_t i = 0; i < t->nbuckets; i++) {
        KmerEntry *e = t->buckets[i];
        while (e) {
            KmerEntry *n = e->next;
            free(e);
            e = n;
        }
    }
    free(t->buckets);
    free(t);
}

static void hash_increment(HashTable *t, const char *ptr) {
    uint64_t h = hash_ptr(ptr, t->k);
    size_t idx = (size_t)(h & (t->nbuckets - 1));
    KmerEntry *e = t->buckets[idx];
    while (e) {
        if (memcmp(e->kmer, ptr, t->k) == 0) {
            e->count++;
            return;
        }
        e = e->next;
    }
    // not found 
    KmerEntry *ne = malloc(sizeof(KmerEntry) + t->k + 1);
    if (!ne) exit(1);
    memcpy(ne->kmer, ptr, t->k);
    ne->kmer[t->k] = '\0';
    ne->count = 1;
    ne->next = t->buckets[idx];
    t->buckets[idx] = ne;
    t->nentries++;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <k>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *input_file = argv[1];
    int k = atoi(argv[2]);
    if (k <= 0) {
        fprintf(stderr, "Error: k must be a positive integer.\n");
        return EXIT_FAILURE;
    }
    k = (k > MAX_KMER) ? MAX_KMER : k;

    FILE *file = fopen(input_file, "rb");
    if (!file) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    if (fseek(file, 0, SEEK_END) != 0) {
        perror("fseek"); fclose(file); return EXIT_FAILURE;
    }
    long file_size = ftell(file);
    if (file_size < 0) { perror("ftell"); fclose(file); return EXIT_FAILURE; }
    if (file_size < k) {
        fclose(file);
        return 0;
    }

    rewind(file);

    char *buffer = malloc((size_t)file_size);
    if (!buffer) { perror("malloc"); fclose(file); return EXIT_FAILURE; }
    size_t read = fread(buffer, 1, (size_t)file_size, file);
    if (read != (size_t)file_size) {
        // proceed with what we read
        file_size = (long)read;
    }
    fclose(file);

    // choose buckets as power of two
    size_t nbuckets = 1u << 16; // 65536 buckets
    HashTable *table = hash_create(k, nbuckets);

    long long startTime = now_ms();
    for (long i = 0; i <= file_size - k; i++) {
        const char *ptr = buffer + i;
        hash_increment(table, ptr);
    }
    long long endTime = now_ms();


    printf("Results:\n");
    for (size_t i = 0; i < table->nbuckets; i++) {
        KmerEntry *e = table->buckets[i];
        while (e) {
            printf("%s: %d\n", e->kmer, e->count);
            e = e->next;
        }
    }
    printf("Time taken: %lld ms\n", endTime - startTime);

    hash_free(table);
    free(buffer);
    return 0;
}