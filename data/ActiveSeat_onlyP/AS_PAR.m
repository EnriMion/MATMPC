function [Par] = AS_PAR()

% %% Input signals assignment
    
    load rif_params.mat
    N = length(rif_accX); %samples number
    Par=[];
    

        for i=1:N
             Par = [Par; rif_accX(i) rif_roll(i) rif_accY(i)];        
        end
   

    %% save your data in the path of your MATMPC
    save([pwd,'\AS_PAR_DATA'], 'Par');
    
    
    clc;
end

