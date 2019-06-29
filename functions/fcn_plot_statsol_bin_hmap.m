function binary_heatmap=fcn_plot_statsol_bin_hmap(stat_sol,prob_thresh,term_verts_inds_cell,nodes,sel_nodes,param_settings,tight_subplot_flag,ranking_flag)

% param_settings=[numsize_plot fontsize hor_gap bottom_marg left_marg];
num_size_plot=param_settings(1); fontsize=param_settings(2); hor_gap=param_settings(3); bottom_marg=param_settings(4); left_marg=param_settings(5);
n=numel(nodes); truth_table_inputs=rem(floor([0:((2^n)-1)].'*pow2(0:-1:-n+1)),2);

% tight_subplot(Nh, Nw, gap, marg_h, marg_w): 
% gap: gaps between the axes in normalized units
% marg_h: margins in height in normalized units
% marg_w: margins in width in normalized units

if ~isempty(prob_thresh)
        term_verts_inds_cell_thresh=arrayfun(@(x) term_verts_inds_cell{x}(stat_sol(term_verts_inds_cell{x})>prob_thresh), 1:numel(term_verts_inds_cell),'un',0);
        % term_verts_inds_cell(stat_sol(term_verts_inds_cell(:))>prob_thresh);
end

if ~isempty(tight_subplot_flag)
    [ha,~]=tight_subplot(numel(term_verts_inds_cell_thresh),1,[hor_gap hor_gap],[bottom_marg 0.01],[left_marg 0.01]);
end

if isempty(sel_nodes); sel_nodes=1:numel(nodes); end

if sum(ismember(arrayfun(@(x) length(term_verts_inds_cell_thresh{x}), 1:numel(term_verts_inds_cell_thresh)),1))==numel(term_verts_inds_cell_thresh)
    [~,ranking]=sort(stat_sol(cell2mat(term_verts_inds_cell_thresh))); term_verts_inds_cell_thresh=term_verts_inds_cell_thresh(flipud(ranking));
end

for k=1:numel(term_verts_inds_cell_thresh)
    
if k==numel(term_verts_inds_cell_thresh); x_ax_leg=nodes(sel_nodes); else x_ax_leg=[]; end
inds=term_verts_inds_cell_thresh{k}; y_ax_leg=round(stat_sol(inds),3); 

    % subplot(numel(term_verts_inds_cell),1,k);
    if ~isempty(tight_subplot_flag)
        axes(ha(k)); 
    else
        subplot(numel(term_verts_inds_cell_thresh),1,k)
    end
    
    % rank by probability
    if ~isempty(ranking_flag)
    [~,ranking]=sort(stat_sol(inds)); ranking=flipud(ranking);
    disp('RANK')
    else
        ranking=1:numel(inds);
    end
    
    binary_heatmap=heatmap(truth_table_inputs(inds(ranking),sel_nodes),x_ax_leg,y_ax_leg(ranking),'%0.0f','TickAngle',90,...
        'Colormap','redblue','MinColorValue',-1,'MaxColorValue',1,'GridLines','-','FontSize',num_size_plot,'ShowAllTicks',true); set(gca,'FontSize',fontsize)
end

hold off