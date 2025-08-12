#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <bitset>
#include <chrono>
#include <omp.h>
#include <io.h>

#define K 13
#define READ_LEN 100
#define MAX_PATH_LEN 260
#define INITIAL_CAPACITY 1000000
#define BASE_TO_BITS(c) (((c) >> 1) & 3)

extern "C"
{
#include "kseq.h"
}

const std::string READ_FOLDER_PATH = "D:\\Data\\Read";
constexpr int KMER_MASK = (1 << (2 * (K - 1))) - 1;
constexpr int BITVECTOR_BIT_LEN = 1 << (2 * K);

static int fread_wrapper(void *fp, void *buf, int len)
{
    return (int)fread(buf, 1, len, (FILE *)fp);
}
KSEQ_INIT(FILE *, fread_wrapper)

class FastqProcessor
{
public:
    std::vector<std::string> reads;
    std::vector<std::string> filepaths;
    double load_time = 0.0;
    double convert_time = 0.0;

    void load_fastq_filepaths(const std::string &folder)
    {
        struct _finddata_t file;
        intptr_t hFile;
        char search_path[MAX_PATH_LEN];
        sprintf_s(search_path, "%s\\*.fastq", folder.c_str());
        hFile = _findfirst(search_path, &file);
        if (hFile == -1)
        {
            std::cerr << "No FASTQ files found." << std::endl;
            exit(1);
        }
        do
        {
            filepaths.emplace_back(folder + "\\" + file.name);
        } while (_findnext(hFile, &file) == 0);
        _findclose(hFile);
    }

    void parallel_load_reads()
    {
        auto start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < filepaths.size(); ++i)
        {
            FILE *fp = fopen(filepaths[i].c_str(), "r");
            if (!fp)
                continue;
            kseq_t *seq = kseq_init(fp);
            while (kseq_read(seq) >= 0)
            {
                if (seq->seq.l < READ_LEN)
                    continue;
                std::string s(seq->seq.s, READ_LEN);
                #pragma omp critical
                reads.push_back(s);
            }
            kseq_destroy(seq);
            fclose(fp);
        }
        auto end = std::chrono::high_resolution_clock::now();
        load_time = std::chrono::duration<double>(end - start).count();
    }

    void convert_to_bitvectors()
    {
        int batch_size;
        if(K <= 7)
            batch_size = 2000000;
        else
            batch_size = 100000;
        for (size_t start_idx = 0; start_idx < reads.size(); start_idx += batch_size)
        {
            size_t end_idx = std::min(start_idx + batch_size, reads.size());
            std::vector<std::bitset<BITVECTOR_BIT_LEN>> batch(end_idx - start_idx);
            auto start = std::chrono::high_resolution_clock::now();
            #pragma omp parallel for schedule(static, 64)
            for (int i = 0; i < static_cast<int>(end_idx - start_idx); ++i)
            {
                const std::string &read = reads[start_idx + i];
                std::bitset<BITVECTOR_BIT_LEN> &bv = batch[i];
                uint32_t val = 0;
                for (int j = 0; j < K; ++j)
                    val = (val << 2) | BASE_TO_BITS(read[j]);
                bv[val] = 1;
                for (int j = 1; j <= READ_LEN - K; ++j)
                {
                    val = ((val & KMER_MASK) << 2) | BASE_TO_BITS(read[j + K - 1]);
                    bv[val] = 1;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            convert_time += std::chrono::duration<double>(end - start).count();
        }
    }

    std::string get_filename(const std::string &path)
    {
        size_t pos = path.find_last_of("\\/");
        return (pos == std::string::npos) ? path : path.substr(pos + 1);
    }

    void report()
    {
        std::cout << " - k = " << K << "\n";
        std::cout << " - Loaded " << reads.size() / 2 << " read pairs from : ";
        if (filepaths.size() >= 2)
            std::cout << get_filename(filepaths[0]) << ", " << get_filename(filepaths[1]) << "\n";
        else if (!filepaths.empty())
            std::cout << get_filename(filepaths[0]) << "\n";

        std::cout << " - FASTQ load time : " << load_time << " s\n";
        std::cout << " - Bitvector convert time : " << convert_time << " s\n";
        std::cout << " - Total runtime : " << (load_time + convert_time) << " s\n";
    }
};

int main()
{
    omp_set_num_threads(omp_get_max_threads());
    FastqProcessor processor;
    processor.load_fastq_filepaths("D:\\Data\\Read");
    processor.parallel_load_reads();
    processor.convert_to_bitvectors();
    processor.report();
    return 0;
}