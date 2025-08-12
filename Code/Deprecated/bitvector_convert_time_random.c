#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <windows.h>
#include <omp.h>

#define READ_LEN 100
#define K 5
#define BITVECTOR_LEN (1 << (2 * K))  // 4^K

unsigned int rand_r_fallback(unsigned int *seed) {
    *seed = (*seed * 1103515245 + 12345) & 0x7fffffff;
    return *seed;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("用法: %s <NUM_READS>\n", argv[0]);
        return 1;
    }

    int NUM_READS = atoi(argv[1]);
    if (NUM_READS <= 0) {
        printf("錯誤：NUM_READS 必須是正整數。\n");
        return 1;
    }

    omp_set_num_threads(8);  // 可以調整為你電腦的核心數
    printf("使用執行緒數: %d\n", omp_get_max_threads());
    printf("開始處理 %d 條序列...\n", NUM_READS);

    LARGE_INTEGER freq, start, end;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);

    #pragma omp parallel for schedule(static)
    for (int r = 0; r < NUM_READS; r++) {
        char sequence[READ_LEN + 1];
        char bitvector[BITVECTOR_LEN] = {0};
        int tid = omp_get_thread_num();
        unsigned int seed = (unsigned int)time(NULL) + r * 7919 + omp_get_thread_num() * 123457;

        for (int i = 0; i < READ_LEN; i++) {
            int base = rand_r_fallback(&seed) % 4;
            sequence[i] = "ACGT"[base];
        }
        // if (r <= 20) {
        //     printf("第 %d\t 條序列：%s\n", r,sequence);
        // }
        sequence[READ_LEN] = '\0';

        int mask = (1 << (2 * K)) - 1;
        int val = 0;

        for (int i = 0; i < K; i++) {
            val <<= 2;
            switch (sequence[i]) {
                case 'A': val += 0; break;
                case 'C': val += 1; break;
                case 'G': val += 2; break;
                case 'T': val += 3; break;
            }
        }
        bitvector[val] = 1;

        for (int i = 1; i <= READ_LEN - K; i++) {
            val <<= 2;
            switch (sequence[i + K - 1]) {
                case 'A': val += 0; break;
                case 'C': val += 1; break;
                case 'G': val += 2; break;
                case 'T': val += 3; break;
            }
            val &= mask;
            bitvector[val] = 1;
        }
    }

    QueryPerformanceCounter(&end);
    double elapsed = (double)(end.QuadPart - start.QuadPart) / freq.QuadPart;

    printf("\n已處理 %d 條序列，每條長度 %d，k = %d\n", NUM_READS, READ_LEN, K);
    printf("每條 bitvector 長度: %d (4^%d)\n", BITVECTOR_LEN, K);
    printf("總轉換耗時: %.9f 秒\n", elapsed);
    printf("平均每條耗時: %.15f 秒\n", elapsed / NUM_READS);

    return 0;
}