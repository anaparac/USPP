function [ W1,W2,H1,H2 ] = PseudoDeflationApproachNMF( X1, X2, k_c, k_d )
    
    alfa = 100;
    beta = 10;
    gama = 10;
    iterations = 100;

    t = 50;     
    k = k_c + k_d;
    valueOfElement = eps;
    
    size1 = size(X1);
    size2 = size(X2);
    
    m = size1(1);
    n1 = size1(2);
    n2 = size2(2);
    
    % inicijalizacija
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
    
    for i = 1 : k
        H1(:,i) = H1(:,i) * norm(W1(:,i));
        H2(:,i) = H2(:,i) * norm(W2(:,i));
        
        W1(:,i) = W1(:,i)/norm(W1(:,i));
        W2(:,i) = W2(:,i)/norm(W2(:,i));
    end
 
    
    R1 = zeros(t, k_d);
    R2 = zeros(t, k_d);
    
    for kds = 0:k_d

       kcs= k_c + k_d - kds;
       for i = 1 : iterations
           
       %--------------------------------------
       % Updateamo Wi-ove
       %--------------------------------------
       
       % update pomoænih podataka  
           H1_ll = H1' * H1;
           H2_ll = H2' * H2;

%            W1c
           for l = 1 : kcs
               if H1_ll(l,l) == 0 
                   H1_ll(l,l) = valueOfElement;
               end
               W1(:, l) =  ( H1_ll(l,l) * W1(:,l) +  X1*H1(:,l) - W1*H1_ll(:,l) + n1  * (alfa / kcs) * W2(:,l) ) / ( H1_ll(l,l) + n1 * alfa );
           end
           W1(W1 < 0) = 0;
           
%            W1d
           if(kds > 0)
            if H1_ll(kds,k) == 0 
               H1_ll(kds,k) = valueOfElement;
            end
            W1(:,k) = W1(:,k) + (X1*H1(:,k) - W1*H1_ll(:,k)- n1*gama/2*W2(:,k))/H1_ll(kds,k);
             if(kds > 1)
                suma = zeros(m,1);
                for p = 1 : (kds - 1)
                    I = zeros(m);
                    for r = 1 : t
                        I(R2(r,p)) = 1;
                    end
                    suma = suma + I*W2(:,k_c+p); 
                end
                
                W1(:,k) = W1(:,k) - n1*beta/(2*(kds-1)) * suma / H1_ll(kds,k);
             end
%             naðemo t najreprezentativnijih kljuènih rijeèi
             [val, ind] = sort(W1(:,k),'descend');
             R2(:,kds) =  sort(ind(1:t));
           end
           
           W1(W1 < 0) = 0;
       % update pomoænih podataka 
           W1_ll = W1' * W1;
           HW1_ll = H1 * W1_ll;
           XW1 = X1' * W1;
     
%            W2c
           for l = 1 : kcs
               if H2_ll(l,l) == 0 
                   H2_ll(l,l) = valueOfElement;
               end              
               W2(:, l) =  ( H2_ll(l,l) * W2(:,l) + X2*H2(:,l) - W2*H2_ll(:,l) + n2  * (alfa / kcs) * W1(:,l) ) / ( H2_ll(l,l) + n2 * alfa );
           end
           
           W2(W2 < 0) = 0;
           
%            W2d
            if(kds > 0)
            if H2_ll(kds,k) == 0 
               H2_ll(kds,k) = valueOfElement;
            end     
            W2(:,k) = W2(:,k) + (X2*H2(:,k) - W2*H2_ll(:,k)- n2*gama/2*W1(:,k))/H2_ll(kds,k);
             if(kds > 1)
                suma = zeros(m,1);
                for p = 1 : (kds - 1)
                    I = zeros(m);
                    for r = 1 : t
                        I(R1(r,p)) = 1;
                    end
                    suma = suma + I*W1(:,k_c+p); 
                end
                
                W2(:,k) = W2(:,k) - n2*beta/(2*(kds-1)) * suma / H2_ll(kds,k);
             end
%             naðemo t najreprezentativnijih kljuènih rijeèi
             [val, ind] = sort(W2(:,k),'descend');
             R1(:,kds) =  sort(ind(1:t));
            end
           
           W2(W2 < 0) = 0;
       % update pomoænih podataka 
           W2_ll = W2' * W2;
           HW2_ll = H2 * W2_ll;
           XW2 = X2' * W2;

           
       %--------------------------------------
       % Updateamo Hi-ove
       %--------------------------------------

           for l =  1 : k
                if W1_ll(l,l) == 0 
                    W1_ll(l,l) = valueOfElement;
                end
               H1(:,l) = H1(:,l) + ( XW1(:,l) - HW1_ll(:,l) ) / W1_ll(l,l);
           end
           
           for l =  1 : k
                if W2_ll(l,l) == 0 
                    W2_ll(l,l) = valueOfElement;
                end               
               H2(:,l) = H2(:,l) + ( XW2(:,l) - HW2_ll(:,l) ) / W2_ll(l,l);
           end

           H1(H1 < 0) = 0;
           H2(H2 < 0) = 0;
           
%           množenje s normom i normiranje           
           for h  = 1 : k  
            pom1 = norm(W1(:,h));
            if norm(W1(:,h)) == 0
                pom1 = eps;
            end      
            H1(:,h) = H1(:,h) * pom1;
            W1(:,h) = W1(:,h) / pom1;
           end
    
           for h  = 1 : k
            pom2 = norm(W2(:,h)); 
            if norm(W2(:,h)) == 0
                pom2 = eps;
            end
            H2(:,h) = H2(:,h) * pom2;
            W2(:,h) = W2(:,h) / pom2;
           end

       end
       
       %izabiremo stupce koje cemo premjestiti na kraj matrica W
       suma = 0;
       N1 = W1(:,1)*H1(:,1)';
       NN1 = N1 - X1;
       NN1(NN1 < 0) = 0;
       mat = N1 - NN1;
       suma  = suma +  norm(mat, 'fro')^2;

       N2 = W2(:,1)*H2(:,1)';
       NN2 = N2 - X2;
       NN2(NN2 < 0) = 0;
       mat = N2 - NN2;
       suma  = suma +  norm(mat, 'fro')^2;
       min_ind = 1;
       min = suma;
       
       for l = 2 : kcs
           suma = 0;
           
           N1 = W1(:,l)*H1(:,l)';
           NN1 = N1 - X1;
           NN1(NN1 < 0) = 0;
           mat = N1 - NN1;
           suma  = suma +  norm(mat, 'fro')^2;
           
           N2 = W2(:,l)*H2(:,l)';
           NN2 = N2 - X2;
           NN2(NN2 < 0) = 0;
           mat = N2 - NN2;
           suma  = suma +  norm(mat, 'fro')^2;
           
           if(suma < min)
               min = suma;
               min_ind = l;
           end          
       end
       
       W1 = [W1(:,1:min_ind - 1), W1(:,min_ind + 1 : k), W1(:,min_ind)];
       H1 = [H1(:,1:min_ind - 1), H1(:,min_ind + 1 : k), H1(:,min_ind)];
       W2 = [W2(:,1:min_ind - 1), W2(:,min_ind + 1 : k), W2(:,min_ind)];
       H2 = [H2(:,1:min_ind - 1), H2(:,min_ind + 1 : k), H2(:,min_ind)];
    end
    
end

