function [PA,EIGENS]=pa_rule_polychoric_missing(MSCORES,missingcode,replicates,smoothing)     

% Code developed by Luis Eduardo Garrido, Ph.D.
% e-mail: garrido.luiseduardo@gmail.com
% This code is based on the methods evaluated and recommended in
% Garrido, L.E., Abad, F.J., & Ponsoda, V. (2013). A new look at Horn's 
% parallel analysis with ordinal variables. Psychological Methods, 18(4), 
% 454-474. 

% INPUT PARAMETERS

% MSCORES: the MATLAB workspace should contain a data matrix named MSCORES 
% with the discrete scores on the variables of interest.

% missingcode: a discrete score used to identify the missing values in the
% MSCORES matrix. Even if there are no missing values a missingcode should
% be given.

% replicates: number of random datasets used to generate the PA criterion
% eigenvalues. 

% smoothing: used to specify if the non-positive definite polychoric
% matrices should be smoothed. A value of 1 indicates that smoothing
% will NOT be conducted, while a value of 2 indicates that smoothing
% will BE conducted. Smoothing is carried out using the Knol & Berger
% eigenvalue method described in Garrido, Abad, and Ponsoda (2013). 

% OUTPUTS

% PA is the number of factors suggested by the method.

% EIGENS is a 2-column matrix that has the real eigenvalues in the first
% column and the criterion eigenvalues in the second column. The criterion
% eigenvalues from the different replicates are combined using the mean 
% statistic. 

% NOTES

% Missing values are handled using pairwise deletion in order to compute 
% the polychoric correlation matrices.  

% The random datasets are generated using random column permutations of the
% real dataset. The missing values are included in the random column
% permutations.

% The syntax file polychoric_proc_missing.m is needed to run this code and
% should be placed in the same working folder. 

% EXAMPLE

% [PA,EIGENS]=pa_rule_polychoric_missing(MSCORES,999,1000,2)   
% This syntax indicates that the matrix of discrete scores is named 
% MSCORES, that the missing values are coded as 999, that 1000
% replicates should be generated, and that the non-positive definite
% matrices should be smoothed. 

% END OF INSTRUCTIONS
                                                    
Size = size(MSCORES);                                                       
N = Size(1,1);                                                                   
Var = Size(1,2);                                                            
                                                                   
% COMPUTING THE REAL EIGENVALUES

MCORR = polychoric_proc_missing(MSCORES,missingcode);                            

if min(eig(MCORR))<=0 && smoothing==2                                       
    [V,D] = eig(MCORR);                                                     
    D = diag(D);                                                            
    D = max(D,0.01);                                                        
    D = diag(D);                                                            
    BB = V*D*V';                                                            
    T = diag(diag(BB))^-0.5;                                                
    MCORR = T*BB*T;                                                           
end
EIGEN = -sort(-eig(MCORR));                                                 

% COMPUTING THE RANDOM EIGENVALUES

Z = zeros(N,Var);                                                           
MEIGENAL = zeros(Var,replicates);                                           
for I=1:replicates
    for J=1:Var                                                             
        X = randperm(N);                                                    
        Y = MSCORES(:,J);                                                   
        Z(:,J) = Y(X);                                                      
    end
    MCORRZ = polychoric_proc_missing(Z,missingcode);
    if min(eig(MCORRZ))<=0 && smoothing==2                                  
        [V,D] = eig(MCORRZ);                                                
        D = diag(D);                                                        
        D = max(D,0.01);                                                    
        D = diag(D);                                                        
        BB = V*D*V';                                                        
        T = diag(diag(BB))^-0.5;                                            
        MCORRZ = T*BB*T;                                                     
    end
    MEIGENAL(:,I) = -sort(-eig(MCORRZ));                                    
    disp(['Replicate = ',num2str(I)])
end
EIGENAL = mean(MEIGENAL,2);                                                 
EIGENS = zeros(Var,2);
EIGENS(:,1) = EIGEN;
EIGENS(:,2) = EIGENAL;

% COMPARING THE REAL AND RANDOM EIGENVALUES

PA = 0;                                                                     
I = 1;                                                                                                                                       
while I<=Var  && EIGEN(I,1)>EIGENAL(I,1)                                     
    PA = PA+1;                                                              
    I = I+1;                                                                
end    
