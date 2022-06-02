function M=PauliM(PauliIndex)
    % Pauli matrices. 
    
    % Translation of a char input. 
    if ischar(PauliIndex)
        switch PauliIndex
            case 'I'
                PauliIndex=0;
            case 'X'
                PauliIndex=1;
            case 'Y'
                PauliIndex=2;
            case 'Z'
                PauliIndex=3;
            otherwise
                error('Invalid PauliIndex!');
        end
    end

    switch PauliIndex
        case 0
            M=speye(2);
        case 1
            M=sparse([0,1;1,0]);
        case 2
            M=sparse([0,-1i;1i,0]);
        case 3
            M=sparse([1,0;0,-1]);
        otherwise
            error('Invalid PauliIndex!');
    end
end