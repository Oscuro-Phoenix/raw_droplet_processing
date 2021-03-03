function [c,ceq] = nlcon(P)
B = [P(1) P(3)/2; 
     P(3)/2 P(2)];  
smalleig = 1e-4; % empirical value to ensure the bilinear form is strictly positive
c = smalleig-eig(B);
ceq = [];
end