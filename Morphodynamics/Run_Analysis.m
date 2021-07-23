%simulate full random walk, pure diffusivity (make all states zero)
%        %Sim_Switching(positions, 11, keep_persist, 1, 0.25, 1);
%simulate balistic (sigma=0) (all states 1)
%        %Sim_Switching(positions, 0, keep_persist, 1, 0, 1);
%simulate pure persistent random walk with different sigma's (all states 1)
%        %Sim_Switching(positions, 0, keep_persist, 1, 0.25, 1);
%simulate switching
%        %Sim_Switching(positions, 5, keep_persist, 1, 0.25, 2);

[shape_positions]=Constrained_Random_Walk(0.1,0.1,2,100,50);
[MSD_shape, MSD_shape_error]=Compute_MSD(shape_positions);
figure(9)
errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_shape, MSD_shape_error);
xlabel('Time',fontsize=20)
ylabel('\sigma^2',fontsize=20)
title('Simulated Shape Space MSD')

D_small = 2;
D_medium = D_small*2;
D_large = D_small*3;
%keep_persist = false;
Pt = 4; %persistence time; bigger = more persistent
Ps = 1.13; %persistence speed; bigger = less persistent?

figure(500)
% random_walk = Sim_Switching(shape_positions, 11, keep_persist, 1, 0.25, 1);
% ballistic_walk = Sim_Switching(shape_positions, 0, keep_persist, 1, 0, 1);
% persist_walk_sig_small = Sim_Switching(shape_positions, 0, keep_persist, 1, sig_small, 1);
% persist_walk_sig_medium = Sim_Switching(shape_positions, 0, keep_persist, 1, sig_medium, 1);
% persist_walk_sig_large = Sim_Switching(shape_positions, 0, keep_persist, 1, sig_large, 1);
% switch_walk = Sim_Switching(shape_positions, 5, keep_persist, 1, sig_small, 1);

random_walk = Sim_Switching_v2(shape_positions, 11, D_small, Pt, Ps);
ballistic_walk = Sim_Switching_v2(shape_positions, 0, 0, Pt, Ps);
persist_walk_sig_small = Sim_Switching_v2(shape_positions, 0, D_small, Pt, Ps);
%persist_walk_sig_medium = Sim_Switching_v2(shape_positions, 0, D_medium, Pt, Ps);
%persist_walk_sig_large = Sim_Switching_v2(shape_positions, 0, D_large, Pt, Ps);
switch_walk = Sim_Switching_v2(shape_positions, 5, D_small, Pt,Ps);

[MSD_random, MSD_error_random] = Compute_MSD(random_walk);
[MSD_ballistic, MSD_error_ballistic] = Compute_MSD(ballistic_walk);
[MSD_persist_s, MSD_error_persist_s] = Compute_MSD(persist_walk_sig_small);
%[MSD_persist_m, MSD_error_persist_m] = Compute_MSD(persist_walk_sig_medium);
%[MSD_persist_l, MSD_error_persist_l] = Compute_MSD(persist_walk_sig_large);
[MSD_switch, MSD_error_switch] = Compute_MSD(switch_walk);

%leginfo={'RW', 'ballistic', sprintf('sig_p=%.2f',D_small), sprintf('sig_p=%.2f',D_medium), sprintf('sig_p=%.2f',D_large), 'switch'};
leginfo={'RW', 'ballistic', 'persist', 'switch'};


errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_random, MSD_error_random);
hold on
errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_ballistic, MSD_error_ballistic);
hold on
errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_persist_s, MSD_error_persist_s);
hold on 
% errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_persist_m, MSD_error_persist_m);
% hold on
% errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_persist_l, MSD_error_persist_l);
% hold on
errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_switch, MSD_error_switch);
xlabel('Time',fontsize=20)
ylabel('\sigma^2',fontsize=20)
title('Simulated Real Space MSD')
legend(leginfo)

figure(600)
errorbar(linspace(0,size(shape_positions,3)-1,size(shape_positions,3))', MSD_switch, MSD_error_switch);
xlabel('Time',fontsize=20)
ylabel('\sigma^2',fontsize=20)
title('Simulated Real Space MSD')

figure(700)
%plot some of the trajectories!
if size(switch_walk,2)==3
    for traj=1:min([5,size(switch_walk,1)])
        xs = squeeze(switch_walk(traj,1,:));
        ys = squeeze(switch_walk(traj,2,:));
        zs = squeeze(switch_walk(traj,3,:));
        plot3(xs,ys,zs,'-x');
        hold on
        
    end
    xlabel('X',fontsize=20)
    ylabel('Y',fontsize=20)
    xlabel('Z', fontsize=20)
elseif size(positions,2)==2
    for traj=1:min([5,size(switch_walk,1)])
        xs = squeeze(switch_walk(traj,1,:));
        ys = squeeze(switch_walk(traj,2,:));
        plot(xs,ys,'-x');
        hold on
    end
    xlabel('X',fontsize=20)
    ylabel('Y',fontsize=20)
else
    disp("Cannot produce plot! Dimensionality too large")
end
title('Simulated Real Space Trajectories')