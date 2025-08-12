// Original algorithm_1 of PIM-quantifier, multi-gene version
Generate index_table for each gene. index_table consists k-mers and associated k_comp classes

Initialize result = ones(m,n) //m is the number of genes, n is the length of k-comp
// This algorithm is in the case of 1 input short read, each short read will call this algorithm to get their own equivalence class
k_mer = short_read[i=0:j=k] // k is the length of k-mer

for k_mer in short_read do  
    for table_len < length(index_table) do
        if k_mer in index_table{table_len} then
            k_comp = index_table{table_len}(k_mer)
            result(table_len,:) = AND(result(table_len,:), k_comp)
        else
            result(table_len,:) = zeros(1,n)
        end if
    end for
    k_mer = short_read[i+1:j+1]
end for
Return result //result indicates the compatible transcripts of all genes