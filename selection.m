 function [ b_cells_trial ] = selection( b_cells_trial, conc, a_act, t_cell_selection) 
% b_cells with affinity higher than a_act can try to internalize the Ag. 
% the probability to internalize the Ag is of a langmuir form. Among these
% successful b_cells, only the best ones are selected by T cells.

%b_cells_trial = b_cells_trial(:, b_cells_trial(:) >= a_act );

for n = 1:size(b_cells_trial,2)
    randoom = rand;
    langmuir = conc*b_cells_trial(n)/(a_act + conc*b_cells_trial(n));
    if randoom >= langmuir
        %b_cell didn't internalize the Ag, so it dies
        b_cells_trial(n) = NaN;
    end  
end
%remove the b_cells that failed to internalize the Ag

b_cells_trial = b_cells_trial(isnan(b_cells_trial) < 1);        

b_cells_trial = b_cells_trial(:, b_cells_trial(:) >= a_act );

%sort the b_cells by increasing value, only the ones at the end are
%selected!
[sorted_b_cells, indexes] = sort(b_cells_trial, 'descend');
n_selected = floor(t_cell_selection*size(b_cells_trial,2));
b_cells_trial = sorted_b_cells(1:n_selected);


%disp(['number of b cells after selection of best performing ' num2str(length(b_cells))]);

end

