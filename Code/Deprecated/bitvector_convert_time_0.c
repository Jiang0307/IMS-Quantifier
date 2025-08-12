#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <omp.h>
#include <io.h>
#include <stdint.h>
#include "kseq.h"

#define K 7
#define READ_LEN 100
#define BITVECTOR_LEN (1 << (2 * K))
#define BITVECTOR_BYTE_LEN (BITVECTOR_LEN / 8)
#define MAX_PATH_LEN 1000
#define BATCH_SIZE 100000
#define KMER_MASK (~(3U << ((K - 1) * 2)))
#define INITIAL_CAPACITY 10000000 // 初始讀入空間，可動態擴增

typedef struct
{
    char **filepaths;
    int count;
} FileList;

// #define KMER_MASK ((1U << (2 * K)) - 1)
static int fread_wrapper(void *fp, void *buf, int len)
{
    return (int)fread(buf, 1, len, (FILE *)fp);
}
KSEQ_INIT(FILE *, fread_wrapper)

int read_count = 0;
int processed = 0;
int batch_id = 0;
double total_elapsed = 0.0;
char **reads = NULL;
// const char *READ_FOLDER_PATH = "C:\\Users\\User\\Desktop\\GitHub\\IMS-quantifier\\Dataset\\Read";
const char *READ_FOLDER_PATH = "D:\\Data\\Read";

// 查表：將 'A', 'C', 'G', 'T' 映射為 00, 01, 10, 11
static const uint8_t char_to_bits[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3};

FileList get_fastq_file_list()
{
    FileList list = {0};
    list.filepaths = malloc(sizeof(char *) * 1000); // 初始容量，後面可改成動態
    list.count = 0;

    struct _finddata_t file;
    char search_path[MAX_PATH_LEN];
    snprintf(search_path, MAX_PATH_LEN, "%s\\*.fastq", READ_FOLDER_PATH);
    intptr_t hFile = _findfirst(search_path, &file);
    if (hFile == -1)
    {
        printf(" - Cannot find any .fastq files\n");
        exit(1);
    }

    do
    {
        char *filepath = malloc(MAX_PATH_LEN);
        snprintf(filepath, MAX_PATH_LEN, "%s\\%s", READ_FOLDER_PATH, file.name);
        list.filepaths[list.count++] = filepath;
    } while (_findnext(hFile, &file) == 0);
    _findclose(hFile);
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
#pragma omp critical // 保證多執行緒不衝突地新增 read
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

int main()
{
    omp_set_num_threads(omp_get_max_threads());
    printf(" - Query performance start\n");

    LARGE_INTEGER freq, t_start, t_end, t_mid;
    QueryPerformanceFrequency(&freq);

    // 測量 FASTQ 載入時間
    printf(" - Load read start\n");
    QueryPerformanceCounter(&t_start);

    FileList files = get_fastq_file_list();
    parallel_load_reads(files);

    QueryPerformanceCounter(&t_mid);
    double load_elapsed = (double)(t_mid.QuadPart - t_start.QuadPart) / freq.QuadPart;
    printf(" - Load FASTQ done, time = %.6f seconds\n", load_elapsed);

    // 測量 bitvector 轉換時間
    // uint8_t *bitvectors = (uint8_t *)_aligned_malloc(BATCH_SIZE * BITVECTOR_BYTE_LEN, 32);
    while (processed < read_count)
    {
        int batch_size = (processed + BATCH_SIZE > read_count) ? (read_count - processed) : BATCH_SIZE;
        size_t bitvector_bytes = (size_t)batch_size * BITVECTOR_BYTE_LEN;
        uint8_t *bitvectors = (uint8_t *)_aligned_malloc(bitvector_bytes, 32); // ✅ 這裡改了
        if (!bitvectors)
        {
            printf(" - Bitvector allocate failed in batch %d\n", batch_id);
            break;
        }
        memset(bitvectors, 0, bitvector_bytes);
        #pragma omp parallel for schedule(guided)
        for (int r = 0; r < batch_size; r++)
        {
            const char *read = reads[processed + r];
            uint8_t *bitvector = bitvectors + r * BITVECTOR_BYTE_LEN;
            uint32_t val = 0;
            for (int i = 0; i < K; i++)
            {
                val = (((val & KMER_MASK) << 2) | char_to_bits[(unsigned char)read[i + K - 1]]) & ((1U << (2 * K)) - 1);
            }
            val &= KMER_MASK;
            _bittestandset((long *)bitvector, val);
            #pragma omp simd
            for (int i = 1; i <= READ_LEN - K; i++)
            {
                val = (((val & KMER_MASK) << 2) | char_to_bits[(unsigned char)read[i + K - 1]]) & ((1U << (2 * K)) - 1);
                _bittestandset((long *)bitvector, val);
            }
        }
        _aligned_free(bitvectors); // ✅ 每個 batch 結束就釋放
        processed += batch_size;
        batch_id++;
    }
    QueryPerformanceCounter(&t_end);
    double convert_elapsed = (double)(t_end.QuadPart - t_mid.QuadPart) / freq.QuadPart;

    // ============================
    // 顯示總時間
    // ============================
    printf(" - Bitvector convert done, time = %.6f seconds\n", convert_elapsed);
    printf(" - Total execution time = %.6f seconds\n", load_elapsed + convert_elapsed);

    return 0;
}