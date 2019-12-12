function [ W1, W2, H1, H2] = BatchProcessingApproachNMF( X1, X2, k_c, k_d )
    
    iterations = 100;
    k = k_c + k_d;
    
    size1 = size(X1);
    size2 = size(X2);
    
    m = size1(1);
    n1 = size1(2);
    n2 = size2(2);
    
    alfa = 100; 
    beta = 10; 
    
    valueOfElement = eps; 
 
    [U1, SIGMA1, V1 ] = svd(X1);
    [U2, SIGMA2, V2 ] = svd(X2);

    W1 =U1(:,1:k)*SIGMA1(1:k,1:k);
    H1 = V1(:,1:k);

    W2 =U2(:,1:k)*SIGMA2(1:k,1:k);
    H2 = V2(:,1:k);

    W2(W2 < 0) = 0;
    W1(W1 < 0) = 0;
    H2(H2 < 0) = 0;
    H1(H1 < 0) = 0;

    for i = 1: k
        if any(W1(:,i) > 0 ) == 0
            W1(1,i) = valueOfElement;
        end
        if any(W2(:,i) > 0 ) == 0
            W2(1,i) = valueOfElement;
        end
        if any(H1(:,i) > 0 ) == 0
            H1(1,i) = valueOfElement;
        end
        if any(H2(:,i) > 0 ) == 0
            H2(1,i) = valueOfElement;
        end
    end

    W1c = W1(:,1:k_c);
    W1d = W1(:, (k_c + 1): k);
    H1c = H1(:,1:k_c);
    H1d = H1(:, (k_c + 1): k);
    W2c = W2(:,1:k_c);
    W2d = W2(:, (k_c + 1): k);
    H2c = H2(:,1:k_c);
    H2d = H2(:, (k_c + 1): k);
    
    
    for i = 1 : k_c
        H1c(:,i) = H1c(:,i) * norm(W1c(:,i));
        H2c(:,i) = H2c(:,i) * norm(W2c(:,i));
        
        W1c(:,i) = W1c(:,i)/norm(W1c(:,i));
        W2c(:,i) = W2c(:,i)/norm(W2c(:,i));
    end
  
    for i = 1 : k_d
        H1d(:,i) = H1d(:,i) * norm(W1d(:,i));
        H2d(:,i) = H2d(:,i) * norm(W2d(:,i));

        W1d(:,i) = W1d(:,i)/norm(W1d(:,i));
        W2d(:,i) = W2d(:,i)/norm(W2d(:,i));
    end
    
   
    % Iteracije   
    for i = 1 : iterations
       W1 = [W1c, W1d];
       W2 = [W2c, W2d];
       H1 = [H1c, H1d];
       H2 = [H2c, H2d];
 
       % update pomocnih varijabli
       H1_ll = H1' * H1;       
       H2_ll = H2' * H2;
      
       %--------------------------------------
       % racunam prvu matricu
       %--------------------------------------
       
       for l = 1 : k_c
           W1c(:,l) =  (H1_ll(l,l) * W1c(:,l) + X1 * H1c(:,l) - W1 * H1_ll(:,l) + n1 * alfa * W2c(:,l) ) /( H1_ll(l,l) + n1 * alfa );     
       end
       W1c(W1c < 0) = 0;
       

       W1(:,1:k_c) = W1c;
       
       
       suma = zeros(m,1);
       for p = 1 : k_d
           suma = suma + W2d(:,p);
       end

       for l = 1 : k_d
        if H1_ll(l+ k_c,l + k_c) == 0 
            H1_ll(l+ k_c,l + k_c) = valueOfElement;
        end
        W1d(:,l) = W1d(:,l) + ( X1 * H1d(:,l ) - W1 * H1_ll(:,l + k_c) - n1 * (beta / 2) * suma)  /  H1_ll(l+ k_c,l + k_c) ;
       end
       W1d(W1d < 0) = 0;
       W1(:,(k_c + 1):k) = W1d;
       
       % update pomocnih varijabli

       W1_ll = W1' * W1;
       HW1_ll = H1 * W1_ll;
       XW1 = X1' * W1; 
       
       for l = 1 : k       
        if W1_ll(l,l) == 0 
            W1_ll(l,l) = valueOfElement;
        end
        H1(:,l) = H1(:,l) + ( XW1(:,l) - HW1_ll(:,l) ) / W1_ll(l,l);            
       end
       H1(H1 < 0) = 0;
       H1c = H1(:,1:k_c);
       H1d = H1(:,(k_c + 1):k);

       
       %--------------------------------------
       % isti postupak za drugu matricu
       %--------------------------------------
            
       for l = 1 : k_c
           W2c(:,l) =  (H2_ll(l,l) * W2c(:,l)+ X2 * H2c(:,l) - W2 * H2_ll(:,l) + n2 * alfa * W1c(:,l))/( H2_ll(l,l) + n2 * alfa );   
       end
       W2c(W2c < 0) = 0;
       W2(:,1:k_c) = W2c;
 
       
       suma = zeros(m,1);
       for p = 1 : k_d
           suma = suma + W1d(:,p);
       end 
       for l = 1 : k_d 
         if H2_ll(l + k_c, l + k_c) == 0 
            H2_ll(l + k_c, l + k_c) = valueOfElement;
        end
       W2d(:,l) = W2d(:,l) + ( X2 * H2d(:,l ) - W2 * H2_ll(:,l + k_c) - n2 * (beta / 2) * suma)  /  H2_ll(l + k_c, l + k_c) ;       
       end
       W2d(W2d < 0) = 0;
       W2(:,(k_c + 1):k) = W2d;

       % update pomocnih varijabli
       
       W2_ll = W2' * W2;
       HW2_ll = H2 * W2_ll;
       XW2 = X2' * W2;      
       
       for l = 1 : k
        if W2_ll(l,l) == 0 
            W2_ll(l,l) = valueOfElement;
        end
         H2(:,l) = H2(:,l) + ( XW2(:,l) - HW2_ll(:,l) ) / W2_ll(l,l);             
       end
       H2(H2 < 0) = 0;
       H2c = H2(:,1:k_c);
       H2d = H2(:,(k_c + 1):k);

       
       % Updateamo vrijednosti : W1c, W1d, W2c, W2d, H1c, H1d, H2c, H2d


        for h  = 1 : k_c
            pom1 = norm(W1c(:,h));
            
            if norm(W1c(:,h)) == 0
                  pom1 = eps;
            end
        
            H1c(:,h) = H1c(:,h) * pom1;
            W1c(:,h) = W1c(:,h)/pom1;
        end
         
        for h  = 1 : k_d
            pom2 = norm(W1d(:,h));
            if norm(W1d(:,h)) == 0
                    pom2 = eps;
            end
        
            H1d(:,h) = H1d(:,h) * pom2;
            W1d(:,h) = W1d(:,h)/pom2;
        end
        
        for h  = 1 : k_c
            pom3 = norm(W2c(:,h));
            if norm(W2c(:,h)) == 0
                  pom3 = eps;
            end
            H2c(:,h) = H2c(:,h) * pom3;
            W2c(:,h) = W2c(:,h)/pom3;
        end
        
        for h  = 1 : k_d
            pom4 = norm(W2d(:,h));
            if norm(W2d(:,h)) == 0
                pom4 = eps;
            end
            H2d(:,h) = H2d(:,h) * pom4;
            W2d(:,h) = W2d(:,h)/pom4;
        end
          
    end 
  
    W1 = [W1c, W1d];
    W2 = [W2c, W2d];
    H1 = [H1c, H1d];
    H2 = [H2c, H2d];
    
   end
    