
xf(:,1)=state_sim(1,:).';
xf2(:,1)=state_sim(1,:).';
for jj=1:2
    if jj==1
        for ii=1:length(output.u)
            sim_input.x = xf(:,ii);
            sim_input.u = output.u(:,ii);
            sim_input.z = input.z(:,ii);
            sim_input.p = input.od(:,ii);

            [xf(:,ii+1), zf] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
            xf(:,ii+1) = full(xf(:,ii+1));
            if ii>=41
                sim_input.x = xf(:,ii+1);
    	        sim_input.u = output.u(:,ii);
                sim_input.z = input.z(:,ii);
                sim_input.p = input.od(:,ii);
        
                [xf(:,ii+1), zf] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
                xf(:,ii+1) = full(xf(:,ii+1));
            end
        end
    elseif jj==2
        for ii=1:length(output.u)
            sim_input.x = xf2(:,ii);
            sim_input.u = output.u(:,ii);
            sim_input.z = input.z(:,ii);
            sim_input.p = input.od(:,ii);

            [xf2(:,ii+1), zf2] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
            xf2(:,ii+1) = full(xf2(:,ii+1));
        end   
    end
end

figure
plot(mem.index_T(:),xf(1,:),'r','Linewidth',1)
hold on
plot(mem.index_T(:),xf2(1,:),'b','Linewidth',1)
title('Position')
legend('Double integration ts>40','Single integration ts>40')
grid on
