function [ res ] = Commonality_norma( W1, W2, k_c)
    matWWc= W1(:, 1:k_c) - W2(:, 1:k_c);
    nfWWc = norm(matWWc, 'fro');
    res = (nfWWc * nfWWc) / k_c; 
end

