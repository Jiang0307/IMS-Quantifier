#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#include <dirent.h>
#include <math.h>
#include "uthash.h"
#include "kseq.h"

#define K 11
#define READ_LEN 100
#define MAX_TRANSCRIPTS 1000
#define MAX_FILES 1024
#define BATCH_SIZE 1000

#define BASE_TO_BITS(c) ((c) == 'A' ? 0 : (c) == 'C' ? 1 : (c) == 'G'   ? 2 : (c) == 'T'   ? 3 : 255)
#define MAX_KMER_VAL (1 << (2 * K))

const char *READ_FOLDER_PATH = "D:\\Data\\Read";
const char *TRANSCRIPT_FOLDER_PATH = "D:\\Data\\Transcript";

static int fread_wrapper(void *fp, void *buf, int len)
{
    return (int)fread(buf, 1, len, (FILE *)fp);
}
KSEQ_INIT(FILE *, fread_wrapper)

typedef struct
{
    char kmer[K + 1];
    int index;
    UT_hash_handle hh;
} KmerEntry;

typedef struct
{
    uint64_t val;
    int index;
    UT_hash_handle hh;
} ValIndex;

KmerEntry *kmer_map = NULL;
int unique_kmer_count = 0;
ValIndex *val_index_map = NULL;
// int *val_to_index = NULL;

char **reads = NULL;
char **transcripts = NULL;
int read_count = 0, transcript_count = 0;

uint8_t **read_bitvectors = NULL;
uint8_t *transcript_bitvectors[MAX_TRANSCRIPTS];
uint8_t base_lookup[256];

// Windows safe strndup
#ifndef _WIN32
#include <string.h>
#else
char *strndup(const char *s, size_t n)
{
    size_t len = strnlen(s, n);
    char *new_str = (char *)malloc(len + 1);
    if (!new_str)
        return NULL;
    memcpy(new_str, s, len);
    new_str[len] = '\0';
    return new_str;
}
#endif

void init_base_lookup()
{
    for (int i = 0; i < 256; i++)
        base_lookup[i] = 255;
    base_lookup['A'] = 0;
    base_lookup['C'] = 1;
    base_lookup['G'] = 2;
    base_lookup['T'] = 3;
}

int is_valid_kmer(const char *kmer)
{
    for (int i = 0; i < K; i++)
    {
        char c = kmer[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
            return 0;
    }
    return 1;
}

void add_kmer_direct(const char *kmer)
{
    if (!is_valid_kmer(kmer))
        return;
    KmerEntry *found;
    HASH_FIND_STR(kmer_map, kmer, found);
    if (!found)
    {
        KmerEntry *entry = malloc(sizeof(KmerEntry));
        strncpy(entry->kmer, kmer, K);
        entry->kmer[K] = '\0';
        entry->index = unique_kmer_count++;
        HASH_ADD_STR(kmer_map, kmer, entry);
    }
}

typedef struct
{
    char *paths[MAX_FILES];
    int count;
} FileList;

void collect_fastq_files(const char *folder, FileList *list)
{
    DIR *dir = opendir(folder);
    struct dirent *entry;
    list->count = 0;
    while ((entry = readdir(dir)) != NULL)
    {
        if (strstr(entry->d_name, ".fastq"))
        {
            list->paths[list->count] = malloc(512);
            snprintf(list->paths[list->count], 512, "%s/%s", folder, entry->d_name);
            list->count++;
        }
    }
    closedir(dir);
}

void load_reads_parallel(const char *folder)
{
    FileList file_list;
    collect_fastq_files(folder, &file_list);
    // printf(" - FASTQ files found : %d\n", file_list.count);

    int capacity = 1000000;
    reads = malloc(capacity * sizeof(char *));

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < file_list.count; i++)
    {
        FILE *fp = fopen(file_list.paths[i], "r");
        if (!fp)
            continue;

        kseq_t *seq = kseq_init(fp);
        while (kseq_read(seq) >= 0)
        {
            if (seq->seq.l < READ_LEN)
                continue;

            for (int j = 0; j <= READ_LEN - K; j++)
            {
                char kmer[K + 1];
                memcpy(kmer, &seq->seq.s[j], K);
                kmer[K] = '\0';
            }

            char *read_copy = strndup(seq->seq.s, READ_LEN);

#pragma omp critical
            {
                if (read_count >= capacity)
                {
                    capacity *= 2;
                    reads = realloc(reads, capacity * sizeof(char *));
                }
                reads[read_count++] = read_copy;
            }
        }
        kseq_destroy(seq);
        fclose(fp);
    }
}

void load_transcripts(const char *folder)
{
    DIR *dir = opendir(folder);
    struct dirent *entry;
    transcripts = malloc(MAX_TRANSCRIPTS * sizeof(char *));

    while ((entry = readdir(dir)) != NULL)
    {
        if (!(strstr(entry->d_name, ".fa") || strstr(entry->d_name, ".fasta")))
            continue;

        char filepath[512];
        snprintf(filepath, sizeof(filepath), "%s/%s", folder, entry->d_name);
        FILE *fp = fopen(filepath, "r");
        if (!fp)
            continue;

        kseq_t *seq = kseq_init(fp);
        while (kseq_read(seq) >= 0 && transcript_count < MAX_TRANSCRIPTS)
        {
            for (int i = 0; i <= seq->seq.l - K; i++)
            {
                char kmer[K + 1];
                memcpy(kmer, &seq->seq.s[i], K);
                kmer[K] = '\0';
            }
            transcripts[transcript_count++] = strdup(seq->seq.s);
        }
        kseq_destroy(seq);
        fclose(fp);
    }
    closedir(dir);
}

void build_kmer_index()
{
    // 掃 read k-mers
    for (int r = 0; r < read_count; r++)
    {
        const char *seq = reads[r];
        for (int i = 0; i <= READ_LEN - K; i++)
        {
            char kmer[K + 1];
            memcpy(kmer, &seq[i], K);
            kmer[K] = '\0';
            add_kmer_direct(kmer);
        }
    }

    // 掃 transcript k-mers
    for (int t = 0; t < transcript_count; t++)
    {
        const char *seq = transcripts[t];
        int len = strlen(seq);
        for (int i = 0; i <= len - K; i++)
        {
            char kmer[K + 1];
            memcpy(kmer, &seq[i], K);
            kmer[K] = '\0';
            add_kmer_direct(kmer);
        }
    }

    // 建 val → index 對應表
    for (KmerEntry *e = kmer_map; e != NULL; e = e->hh.next)
    {
        uint32_t val = 0;
        for (int i = 0; i < K; i++)
        {
            val = (val << 2) | base_lookup[(unsigned char)e->kmer[i]];
        }
        ValIndex *v = malloc(sizeof(ValIndex));
        v->val = val;
        v->index = e->index;
        HASH_ADD(hh, val_index_map, val, sizeof(uint64_t), v);
    }
}

void build_read_bitvectors()
{
    int bytes = (unique_kmer_count + 7) / 8;
    read_bitvectors = malloc(sizeof(uint8_t *) * read_count);
    uint8_t *bitvector_pool = calloc(read_count * bytes, 1);

#pragma omp parallel for
    for (int r = 0; r < read_count; r++)
    {
        const char *seq = reads[r];
        uint8_t *bv = bitvector_pool + r * bytes;
        read_bitvectors[r] = bv;

        uint32_t val = 0;
        int valid = 1;

        // encode first K
        for (int i = 0; i < K; i++)
        {
            uint8_t b = base_lookup[(unsigned char)seq[i]];
            if (b == 255)
            {
                valid = 0;
                break;
            }
            val = (val << 2) | b;
        }

        if (valid)
        {
            ValIndex *found;
            HASH_FIND(hh, val_index_map, &val, sizeof(uint64_t), found);
            if (found)
            {
                int idx = found->index;
                bv[idx / 8] |= 1 << (idx % 8);
            }

            uint32_t mask = (1 << (2 * K)) - 1;
            for (int i = K; i < READ_LEN; i++)
            {
                uint8_t b = base_lookup[(unsigned char)seq[i]];
                if (b == 255)
                {
                    break;
                } // fail fast
                val = ((val << 2) | b) & mask;

                HASH_FIND(hh, val_index_map, &val, sizeof(uint64_t), found);
                if (found)
                {
                    int idx = found->index;
                    bv[idx / 8] |= 1 << (idx % 8);
                }
            }
        }
    }

    printf(" - Read bitvectors built (%d bytes each)\n", bytes);
}

void build_transcript_bitvectors()
{
    int bytes = (unique_kmer_count + 7) / 8;

#pragma omp parallel for
    for (int t = 0; t < transcript_count; t++)
    {
        const char *seq = transcripts[t];
        int len = strlen(seq);
        uint8_t *bv = calloc(bytes, 1);
        for (int i = 0; i <= len - K; i++)
        {
            uint32_t val = 0;
            int valid = 1;
            for (int j = 0; j < K; j++)
            {
                uint8_t bits = BASE_TO_BITS(seq[i + j]);
                if (bits == 255)
                {
                    valid = 0;
                    break;
                }
                val = (val << 2) | bits;
            }
            if (valid)
            {
                ValIndex *found;
                HASH_FIND(hh, val_index_map, &val, sizeof(uint64_t), found);
                if (found)
                {
                    int idx = found->index;
                    bv[idx / 8] |= 1 << (idx % 8);
                }
            }
        }
        transcript_bitvectors[t] = bv;
    }
    // printf(" - Transcript bitvectors built (%d bytes each)\n", bytes);
}

int main()
{
    omp_set_num_threads(omp_get_max_threads());
    init_base_lookup();

    load_transcripts(TRANSCRIPT_FOLDER_PATH);
    double load_start = omp_get_wtime();
    load_reads_parallel(READ_FOLDER_PATH);
    double load_end = omp_get_wtime();
    printf(" - Loaded %d reads, %d transcripts\n", read_count, transcript_count);
    printf(" - Load time : %.4f seconds\n", load_end - load_start);

    double convert_start = omp_get_wtime();
    build_kmer_index();
    build_read_bitvectors();
    double convert_end = omp_get_wtime();
    build_transcript_bitvectors();

    printf("\n - k = %d\n", K);
    printf(" - Bitvector length before compaction : %.0f\n", pow(4, K));
    printf(" - Bitvector length after compaction  : %d\n", unique_kmer_count);
    printf(" - Convert time : %.4f seconds\n", convert_end-convert_start);
    printf(" - Total time : %.4f seconds\n", (load_end-load_start) + (convert_end-convert_start));

    return 0;
}