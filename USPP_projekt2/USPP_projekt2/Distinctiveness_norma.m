function [ sum ] = Distinctiveness_norma( W1, W2, k_d)
    W1(W1 == 0) = eps;
    W2(W2 == 0) = eps;
    sum = 0;
    for  i = 1 : k_d
        for j = 1 : k_d
    
            log1 = log10(W1(:, (k_c + i))); 
            log2 = log10(W2(:, (k_c + j)));
                      
            term = W1(:, (k_c + i))' * log1 ...
                 + W2(:, (k_c + j))' * log2 ...
                 - W1(:, (k_c  + i))' * log2 ...
                 - W2(:, (k_c + j))' * log1;
          
            sum = sum + term;
                     
        end
    end
    sum = sum / (k_d * k_d * 2);
end


