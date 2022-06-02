function M=OnsitePauliMString(PauliIndices,SiteIndex,N_Sites)
    %----------------------------------------------------------------------
    % Construct a string of Pauli matrices specified by PauliIndices acting
    % on a given site specified by the SiteIndex. N_Sites is the total 
    % number of sites. 
    %----------------------------------------------------------------------
    StringLength=length(PauliIndices);
    
    assert((1<=SiteIndex)&&(SiteIndex<=N_Sites),'Invalid SiteIndex!');
    assert(StringLength<=N_Sites,'The Pauli matrix string is too long!');
    
    if (SiteIndex+StringLength-1)<=N_Sites
        M=kron(kron(speye(2^(SiteIndex-1)),PauliMString(PauliIndices)),speye(2^(N_Sites-(SiteIndex+StringLength-1))));
    else
        % Need to break into two parts. 
        LengthFirstPart=N_Sites-SiteIndex+1; % Length of the first part. 
        
        M=OnsitePauliMString(PauliIndices(1:LengthFirstPart),SiteIndex,N_Sites)*...
          OnsitePauliMString(PauliIndices((LengthFirstPart+1):StringLength),1,N_Sites);
    end
end