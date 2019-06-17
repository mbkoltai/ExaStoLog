% function that plots A (trans matrix), K (kinetic matrix) and solutions
function fcn_plot_A_K_stat_sol(A, nodes, sel_nodes, stat_sol, x0, min_max_col,fontsize,barwidth_states_val,nonzero_flag)

[stationary_node_vals,init_node_vals]=fcn_calc_init_stat_nodevals(x0,stat_sol);
min_col=min_max_col(1); max_col=min_max_col(2);
if isempty(sel_nodes)
    sel_nodes=1:numel(nodes); 
end

if ~isempty(A)

fontsize_hm=fontsize(1); fontsize_stat_sol=fontsize(2);
% min_col=-1; max_col=1;fontsize_hm=10;
% transition matrix
fig_subpl1=subplot(1,3,1); set(fig_subpl1,'Position',[0.032 0.11 0.32 0.82]); % [0.032 0.11 0.32 0.82]
if (issparse(A))
    spy(A)
else
    heatmap(A,1:2^length(nodes),1:2^length(nodes),'%0.2f','TickAngle',90,'MinColorValue',min_col,'MaxColorValue',max_col,'Colormap','redblue',... 
        'GridLines', '-', 'FontSize', fontsize_hm, 'ShowAllTicks', true); 
end
% title 
if full(min(A(:)))==0
    title('A (transition matrix)', 'FontSize', 20, 'FontWeight','normal');
else
    title('K (kinetic matrix)', 'FontSize', 20, 'FontWeight','normal');
end

subplot(1,3,2); c=strsplit(num2str(1:2^length(nodes))); 
nnz_vals=[find(stat_sol>0); find(x0>0.05)]; 
plot_vals=[x0(nnz_vals)'; stat_sol(nnz_vals)']';
bar_subpl2=barh(1:numel(nnz_vals),plot_vals,'grouped','BarWidth',1); xlim([0 1.02]); set(gca,'yticklabel',num2cell(nnz_vals))
set(bar_subpl2(size(plot_vals,2)),'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'BarWidth',barwidth_states_val); 
% ylim([0 round(numel(x0)*1.02)]); 
grid on; 
title('stationary solution: states','FontSize',fontsize_stat_sol,'FontWeight','normal'); % set(gca,'YtickLabel',{c{:},'sum'}); 
legend({'x_0', 'steady state'},'Location','SouthEast','FontSize',fontsize_stat_sol*1.2);

subplot(1,3,3); bar_subpl3=barh(flipud([init_node_vals(sel_nodes); stationary_node_vals(sel_nodes)]'), 'grouped'); grid on; 
title('stationary solution: nodes','FontSize',fontsize_stat_sol,'FontWeight','normal'); 
set(gca,'YtickLabel',strrep(fliplr(nodes(sel_nodes)),'_',' ')); shift=0.4; ylim([1-shift numel(sel_nodes)+shift])
set(bar_subpl3(2),'FaceColor',[1 0 0],'EdgeColor',[1 0 0]);

else % without the matrix
    
if ~isempty(nonzero_flag)
    nnz_vals=find(stat_sol>0); state_vals=stat_sol(nnz_vals);
else
    state_vals=stat_sol;
end

fig_subpl1=subplot(1,2,1); barh(state_vals,'BarWidth',barwidth_states_val,'FaceColor',[0 0.5 0],'EdgeColor',[0 0.5 0]); 
disp_lim=0.02; xlim([0 1]);  % set(fig_subpl1,'Position',[0.06 0.11 0.4 0.8]); 
if ~isempty(nonzero_flag); set(gca,'yticklabel',num2cell(nnz_vals)); else
set(fig_subpl1,'ytick','') % ylim([find(state_vals>disp_lim,1,'first')-0.05*numel(state_vals) find(state_vals>disp_lim,1,'last')+0.05*numel(state_vals)]); 
end

fontsize=fontsize(1);
title('states','FontWeight','normal'); xlabel('stationary probability'); set(gca, 'FontSize', fontsize)
fig_subpl1.XGrid = 'on'; % % set(fig_subpl1,'ytick',find(stat_sol>disp_lim)); 

fig_subpl2=subplot(1,2,2); % node_names = {'dna\_dam','chek1','mk2','atm\_atr','hr','cdc25b','g2m','cell\_death'}; 
if ~isempty(init_node_vals)
    plot_vals=flipud([init_node_vals(sel_nodes); stationary_node_vals(sel_nodes)]');
else
    plot_vals=fliplr(stationary_node_vals(sel_nodes))';
end
bar1=barh(plot_vals,'grouped'); set(gca, 'YTickLabel',fliplr(strrep(nodes(sel_nodes),'_',' '))); 
set(gca,'FontSize',fontsize); xlabel('stationary probability'); n=size(plot_vals,2); set(bar1(n),'FaceColor',[1 0 0]);
xlim([0 1.05]); set(fig_subpl2,'Position',[0.56 0.11 0.4 0.8]); ylim([0.5 numel(sel_nodes)+0.5]); grid on
if n>1
legend({'x_0', 'steady state'},'Location','SouthEast','FontSize',fontsize);
end
    
end