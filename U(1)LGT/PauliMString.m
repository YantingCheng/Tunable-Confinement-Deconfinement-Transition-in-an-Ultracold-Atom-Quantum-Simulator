function M=PauliMString(PauliIndices)
    StringLength=length(PauliIndices);
    M=1;
    for i=1:StringLength
        M=kron(M,PauliM(PauliIndices(i)));
    end
end