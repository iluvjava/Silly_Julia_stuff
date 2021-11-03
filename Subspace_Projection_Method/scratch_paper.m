N = 5;
T1 = rand(N); 
T1 = T1*T1'; 
T1 = T1 - triu(T1, 2) - tril(T1, -2);
%% LDL: 
[L, D] = ldl(T1);


