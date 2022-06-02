function H=SpinChainH(PauliIndices,N_Sites,LeadingSiteList)
    %----------------------------------------------------------------------
    % Construct a translation invariant spin Hamiltonian. 
    % PauliIndices is a vector specifying the terms to construct. E.g. if
    % PauliIndices == [3,1,3], then the terms Z_i X_{i+1} Z_{i+2} will be
    % contructed. 
    % N_Sites is the total number of sites. 
    % 
    % LeadingSiteList is a vector specifying a set of sites. For each site
    % index i in this list, there will be a term whose leading index is i. 
    % In the example above, if LeadingSiteList=2:2:10, the Hamiltonian will
    % contain terms Z_2 X_3 Z_4 + Z_4 X_5 Z_6 ... + Z_10 X_11 Z_12. 
    % 
    % LeadingSiteList can also be a char string specifying the boundary
    % condition: 'obc' or 'pbc'. 
    %----------------------------------------------------------------------
    TermLength=length(PauliIndices);
    Dim=2^N_Sites; % Hilbert space dimension. 
    
    if nargin<3
        LeadingSiteList='obc';
    end
    
    if ischar(LeadingSiteList)
        if strcmp(LeadingSiteList,'obc')
            LeadingSiteList=1:(N_Sites-TermLength+1);
        elseif strcmp(LeadingSiteList,'pbc')
            LeadingSiteList=1:N_Sites;
        else
            error('Invalid boundary condition!');
        end
    end    
    
    H=sparse(Dim,Dim); 
    for i_Site=LeadingSiteList
        H=H+OnsitePauliMString(PauliIndices,i_Site,N_Sites);
    end
end