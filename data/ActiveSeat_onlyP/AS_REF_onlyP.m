function [Yref] = AS_REF_onlyP(flag_hp)

    %% find your path to the original active seat model files

% %% Input signals assignment
    
    load rif_pressione.mat
    
if flag_hp == 0
    N = length(rif_pressione);
    Yref=[];


    for i=1:N
            Yref = [Yref; rif_pressione(i), 0 ];
    end

end

if flag_hp == 1
    N = length(rif_pressione_hp);
    Yref=[];


    for i=1:N
            Yref = [Yref; rif_pressione_hp(i), 0 ];
    end

end

    %% save your data in the path of your MATMPC
    save([pwd,'\AS_REF_DATA_onlyP'], 'rif_pressione');    
    
    clc;
end

