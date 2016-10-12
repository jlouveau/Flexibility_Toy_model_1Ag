function [ pop_time, total_exit_cells, neutralized, breadth ] = analysis( number_recycled_b_cells, number_exit_cells, exit_cells, n_trial_max, a_act, n_cycle_max, p_mut, p_recycle, t_cell_selection)
%number_recycled_b_cells, exit_cells, n_trial_max, a_act, n_cycle_max, p_mut, p_recycle, t_cell_selection )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% population of GC B cells
pop_time = mean(number_recycled_b_cells,1);

figure(); plot(pop_time);
title({['Population of GC b cells over time for a single antigen with mutations only in the CDR']; ['averaged over ', num2str(n_trial_max), ' trials']; [' with proba recycle = ', num2str(p_recycle) ', proba mutation = ' num2str(p_mut) ' and t cell selection rate = ' num2str(t_cell_selection)]});
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6)
figure(); 
for i = 1:n_trial_max
    plot(number_recycled_b_cells(i,:)); hold on;
end
title({['Population of GC b cells over time for a single antigen with mutations only in the CDR']; [' with proba recycle = ', num2str(p_recycle) ', proba mutation = ' num2str(p_mut) ' and t cell selection rate = ' num2str(t_cell_selection)]});
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6)
%% breadth 
% fraction of exit_cells with affinity for both Ags above a threshold for different thresholds
% exit_cells = zeros(n_trial_max, n_cycle_max, floor(n_max_Bcells/4));

total_exit_cells = mean(number_exit_cells,1);
thresholds = linspace(a_act, 5+a_act,5);
neutralized = zeros(n_trial_max, n_cycle_max,length(thresholds));
breadth = zeros(n_cycle_max,length(thresholds));

for t = 1:length(thresholds)
    for i = 1:size(exit_cells,1)
        for j = 1:size(exit_cells,2)
            for k = 1:size(exit_cells,3)
                   if exit_cells(i,j,k) > thresholds(t)
                        neutralized(i,j,t) = neutralized(i,j,t) +1;
                   end
            end
        end
    end
end

neutralized_mean = mean(neutralized,1);

figure();
plot(total_exit_cells);
title('Number exit cells');
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6)
%title({['Number of b cells that exit the GC over time for a single antigen with mutations only in the CDR'];  ['averaged over ', num2str(n_trial_max), ' trials']; [' with proba recycle = ', num2str(p_recycle) ', proba mutation = ' num2str(p_mut) ' and t cell selection rate = ' num2str(t_cell_selection)]});

figure();
for t = 1:length(thresholds)
    plot(neutralized_mean(:,:,t)); hold on;
end
title('Neutralized');
legendCell = strcat(strtrim(cellstr(num2str(thresholds(:)))));
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6)
legend(legendCell,'fontsize',6, 'Position', [0.75,0.65,0.25,0.25]);

for c = 1:n_cycle_max
    if total_exit_cells(c) == 0
        breadth(c,t) = 0;
    else
        for t = 1:length(thresholds) 
            breadth(c,t) = neutralized_mean(1,c,t);
        end
        breadth(c,:) = breadth(c,:)/total_exit_cells(c);
    end
end

figure();
for t = 1:length(thresholds)
    plot(breadth(:,t)); hold on;
end
title('Breadth');
xlabel('Number of cycles', 'Fontweight', 'bold');
set(gca,'FontSize',6)
legend(legendCell, 'fontsize',6, 'Position', [0.65,0.65,0.15,0.25]);


end
