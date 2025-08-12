#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <dirent.h>
#include <sys/time.h>
#include <unistd.h>

#include "kseq.h"

#define K 15
#define BASE_TO_BITS(c) (((c) >> 1) & 3)
#define KMER_MASK ((1U << (2 * (K - 1))) - 1)
#define BITVECTOR_BIT_LEN (1U << (2 * K))
#define BITVECTOR_BYTE_LEN (BITVECTOR_BIT_LEN / 8)
#define READ_LEN 100
#define MAX_PATH_LEN 300
#define INITIAL_CAPACITY 1000000
#define BATCH_SIZE 5000

const char *READ_FOLDER_PATH = "D:\\Data\\Read";

char **reads = NULL;
int read_count = 0;
double load_time = 0.0;
double convert_time = 0.0;

static int fread_wrapper(void *fp, void *buf, int len)
{
    return (int)fread(buf, 1, len, (FILE *)fp);
}
KSEQ_INIT(FILE *, fread_wrapper)

static const uint8_t char_to_bits[128] = 
{
    ['A'] = 0,
    ['C'] = 1,
    ['G'] = 2,
    ['T'] = 3,
    ['a'] = 0,
    ['c'] = 1,
    ['g'] = 2,
    ['t'] = 3,
};

typedef struct
{
    char **filepaths;
    int count;
} FileList;

double get_time_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

FileList get_fastq_file_list()
{
    FileList list = {0};
    list.filepaths = malloc(sizeof(char *) * 1000);
    list.count = 0;

    DIR *dir = opendir(READ_FOLDER_PATH);
    if (!dir)
    {
        perror("opendir");
        exit(1);
    }

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL)
    {
        if (strstr(entry->d_name, ".fastq"))
        {
            char *filepath = malloc(MAX_PATH_LEN);
            snprintf(filepath, MAX_PATH_LEN, "%s/%s", READ_FOLDER_PATH, entry->d_name);
            list.filepaths[list.count++] = filepath;
        }
    }
    closedir(dir);
    return list;
}

void parallel_load_reads(FileList file_list)
{
    int capacity = INITIAL_CAPACITY;
    reads = malloc(sizeof(char *) * capacity);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < file_list.count; i++)
    {
        FILE *fp = fopen(file_list.filepaths[i], "r");
        if (!fp)
            continue;

        kseq_t *seq = kseq_init(fp);
        while (kseq_read(seq) >= 0)
        {
            if (seq->seq.l < READ_LEN)
                continue;
            #pragma omp critical
            {
                if (read_count >= capacity)
                {
                    capacity *= 2;
                    reads = realloc(reads, sizeof(char *) * capacity);
                }
                reads[read_count] = malloc(READ_LEN + 1);
                strncpy(reads[read_count], seq->seq.s, READ_LEN);
                reads[read_count][READ_LEN] = '\0';
                read_count++;
            }
        }
        kseq_destroy(seq);
        fclose(fp);
    }
}

void convert_reads_to_bitvectors_in_batches()
{
    static uint8_t *bitvectors = NULL;
    int num_batch = BATCH_SIZE;
    if (!bitvectors)
    {
        if(K >= 15)
            num_batch = BATCH_SIZE / 40;
        bitvectors = (uint8_t *)calloc(num_batch, BITVECTOR_BYTE_LEN);
        if (!bitvectors)
        {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
    }

    for (int start = 0; start < read_count; start += num_batch)
    {
        int end = (start + num_batch > read_count) ? read_count : start + num_batch;
        int batch_size = end - start;

        double convert_start = get_time_sec();

        if (batch_size < num_batch)
        {
            memset(bitvectors, 0, num_batch * BITVECTOR_BYTE_LEN);
        }

        #pragma omp parallel for schedule(static)
        for (int r = 0; r < batch_size; ++r)
        {
            const char *read = reads[start + r];
            uint8_t *bitvector = bitvectors + r * BITVECTOR_BYTE_LEN;

            uint32_t val = 0;
            for (int i = 0; i < K; ++i)
            {
                val = (val << 2) | BASE_TO_BITS(read[i]);
            }
            bitvector[val >> 3] |= 1 << (val & 7);

            for (int i = 1; i <= READ_LEN - K; ++i)
            {
                val = ((val & KMER_MASK) << 2) | BASE_TO_BITS(read[i + K - 1]);
                bitvector[val >> 3] |= 1 << (val & 7);
            }
        }

        double convert_end = get_time_sec();
        convert_time += convert_end - convert_start;
    }
}

const char *get_filename(const char *path)
{
    const char *slash = strrchr(path, '/');
    return slash ? slash + 1 : path;
}

int main()
{
    omp_set_num_threads(omp_get_max_threads());

    double load_start = get_time_sec();
    FileList file_list = get_fastq_file_list();
    parallel_load_reads(file_list);
    double load_end = get_time_sec();
    load_time = load_end - load_start;

    convert_reads_to_bitvectors_in_batches();

    printf(" - k = %d\n", K);
    if (file_list.count >= 2)
        printf(" - Loaded %d read pairs from %s, %s\n", read_count / 2, get_filename(file_list.filepaths[0]), get_filename(file_list.filepaths[1]));
    else
        printf(" - Loaded %d reads from %s\n", read_count, file_list.count > 0 ? get_filename(file_list.filepaths[0]) : "no files");
    printf(" - FASTQ load time : %.4f s\n", load_time);
    printf(" - Bitvector convert time : %.4f s\n", convert_time);
    printf(" - Total runtime : %.4f s\n", (load_time + convert_time));

    free(reads);
    return 0;
}