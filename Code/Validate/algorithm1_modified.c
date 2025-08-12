// Modified algorithm_1 of PIM-quantifier, single-gene version
Generate index_table for 1 gene. Index_table consists k-mers and associated k_comp classes

Initialize result = ones(n) // n is the #transcript
k_mer = short_read[i=0:j=k] // k is the length of k-mer

for k_mer in short_read do
    if k_mer in index_table then
        k_comp = index_table(k_mer)
        result = AND(result, k_comp)
    else
        result = zeros(n)
    end if
    k_mer = short_read[i+1:j+1]
end for
return result // result indicates the equivalence class