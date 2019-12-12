function [ res ] = Reconstruciton_error( X1, X2, W1, W2, H1, H2)
    X1_ = W1*H1';
    X2_ = W2*H2';
    mat1 = X1- X1_;
    mat2 = X2- X2_;

    size1 = size(X1);
    n1 = size1(2);

    size2 = size(X2);
    n2 = size2(2);

    nf1 =  norm(mat1, 'fro');
    nf2 =  norm(mat2, 'fro');

    res = 1/n1 * nf1 * nf1 + 1/n2 * nf2 * nf2; 
end
