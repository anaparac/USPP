function [ X1, X2, W1, W2, H1, H2 ] = synteticMatrix( )

common = 6; 
discriminative = 4; 
k = common + discriminative; 

m = 1400;
n = 300; 

W1C = zeros(m, common);
W2C = zeros(m, common);
W1D = zeros(m, discriminative);
W2D = zeros(m, discriminative);

H1C = zeros(n, common);
H2C = zeros(n, common);
H1D = zeros(n, discriminative);
H2D = zeros(n, discriminative);

for l= 1:common
    for p=1:m
        if p > 100*(l-1) && p <= 100*l
            W1C(p,l) = 1;
            W2C(p,l) = 1;
        end
    end
end


for l= 1:discriminative
    for p=1:m
        idx = 600 + 100*(l-1);
        if p > idx && p <= idx + 100;
            W1D(p,l) = 1;
        end
    end
end

for l= 1:discriminative
    for p=1:m
        idx = 1000 + 100*(l-1);
        if p > idx && p <= idx + 100;
            W2D(p,l) = 1;
        end
    end
end

for i = 1 : n
    stup1 = randi(k); 
    if stup1 <= common
        H1C(i, stup1) = 1;
    else 
        stup1 = stup1 - common; 
        H1D(i, stup1) = 1; 
    end
    
    stup2 = randi(k);
    if stup2 <= common
        H2C(i, stup2) = 1;
    else 
        stup2 = stup2 - common; 
        H2D(i, stup2) = 1; 
    end
end


W1 = [W1C, W1D]; 
W2 = [W2C, W2D];
H1 = [H1C, H1D];
H2 = [H2C, H2D];

awgnRand = 90;%randi([80,100]);  
W1 = awgn ( W1, awgnRand);

awgnRand = 90;%randi([80,100]);
W2 = awgn ( W2, awgnRand); 

awgnRand = 90;%randi([80,100]);
H1 = awgn ( H1, awgnRand); 

awgnRand = 90;%randi([80,100]);
H2 = awgn ( H2, awgnRand);


W1(W1<0) = abs(W1(W1<0));
W2(W2<0) = abs(W2(W2<0));

H1(H1<0) = abs(H1(H1<0));
H2(H2<0) = abs(H2(H2<0));

X1 = W1 * H1'; 
X2 = W2 * H2';

awgnRand = randi([80,100]);
X1 = awgn (X1, awgnRand);

awgnRand = 80;%randi([80,100]);
X2 = awgn (X2, awgnRand);


X1(X1<0) = abs(X1(X1<0));
X2(X2<0) = abs(X2(X2<0));
end

