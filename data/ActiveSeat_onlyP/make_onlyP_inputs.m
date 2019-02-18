%mi crea il vettore di parametri (provenienti dal vestibolare) dal file MPC sul sistema globale

function make_onlyP_inputs()
    load data_MPC_ActiveSeat_full.mat
    rif_daccX =[];
    rif_daccY = [];
    rif_accX =[0];
    rif_accY =[0];
    rif_roll = [0];

    for i=1: length(controls_MPC)
        rif_daccX = [rif_daccX, controls_MPC(i,1)];
        rif_daccY = [rif_daccY, controls_MPC(i,3)];
    
    end

    for i=2: length(state_sim)
        rif_accX = [rif_accX,state_sim(i,8)];
        rif_accY = [rif_accY,state_sim(i,17)];
        rif_roll = [rif_roll,state_sim(i,14)];
    
    end

    save([pwd,'\rif_input_pres'],'rif_accX','rif_accY','rif_roll','rif_daccX','rif_daccY');
end
