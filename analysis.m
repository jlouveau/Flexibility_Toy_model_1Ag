function [ pop_time ] = analysis( number_recycled_b_cells, n_trial_max, p_mut, p_recycle, t_cell_selection )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

pop_time = mean(number_recycled_b_cells,1);

figure(); plot(pop_time);
title({['Population of GC b cells over time for a single antigen with mutations only in the CDR'], ['averaged over ', num2str(n_trial_max), ' trials']; [' with proba recycle = ', num2str(p_recycle) ', proba mutation = ' num2str(p_mut) ' and t cell selection rate = ' num2str(t_cell_selection)]});

figure(); 
for i = 1:n_trial_max
    plot(number_recycled_b_cells(i,:)); hold on;
end
title({['Population of GC b cells over time for a single antigen with mutations only in the CDR']; [' with proba recycle = ', num2str(p_recycle) ', proba mutation = ' num2str(p_mut) ' and t cell selection rate = ' num2str(t_cell_selection)]});
end

