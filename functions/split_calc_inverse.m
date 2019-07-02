function [stat_sol_blocks,term_verts_cell,cell_subgraphs]=split_calc_inverse(A_sparse,transition_rates_table,x0)

% is the STG disconnected?
stat_sol_blocks=zeros(numel(x0),1);
A_digraph = digraph(A_sparse,'omitselfloops'); subnetws=conncomp(A_digraph,'Type','weak'); num_subnets=length(unique(subnetws));
% preallocate cell of term vertices and of subgraphs
term_verts_cell = cell(num_subnets,1); cell_subgraphs=cell(num_subnets,1);

for counter=1:num_subnets

submatrix_inds=find(subnetws==counter);
cell_subgraphs{counter}=submatrix_inds;

if sum(x0(submatrix_inds))>0

A_sparse_sub=A_sparse(submatrix_inds,submatrix_inds);
dim_matr=size(A_sparse_sub,1); scc_submat=conncomp(digraph(A_sparse_sub,'omitselfloops'),'Type','strong');

% IF all SCCs are single vertices
if numel(unique(scc_submat))==dim_matr
    % disp('no (transient or terminal) cycles')

% function to reorder vertices and keep ordering
[A_sparse_sub_reordered_terminal,sorted_vertices_terminal_bottom]=fcn_toposort_terminals_bottom(A_sparse_sub);
K_sp_sub_reord = (A_sparse_sub_reordered_terminal' - speye(dim_matr,dim_matr))*sum(transition_rates_table(:));

stat_sol_submatr_blocks=fcn_block_inversion(K_sp_sub_reord,sorted_vertices_terminal_bottom,x0,submatrix_inds);
% global solution: stat_sol_split(rel_states)=stat_sol_submatr;
stat_sol_blocks(submatrix_inds(sorted_vertices_terminal_bottom))=stat_sol_submatr_blocks;
term_verts_cell{counter} = num2cell(intersect(submatrix_inds,find(stat_sol_blocks>0))); % find(stat_sol_blocks>0)
% disp(counter); disp(find(stat_sol_submatr_blocks>0)'); 
% disp(intersect(submatrix_inds,find(stat_sol_blocks>0))') % reshape(sort(submatrix_inds(sorted_vertices_terminal_bottom)),128,64)

else % if there are cycles we need reordering of metagraph of SCC
    disp('cycles in STG');
    
    % if entire graph is one connected component, no reordering needed
    if numel(unique(scc_submat))==1
        K_sp_sub_reord = (A_sparse_sub' - speye(dim_matr,dim_matr))*sum(transition_rates_table(:));
        kernel_col=((-1)^(dim_matr-1))*fcn_adjug_matrix(K_sp_sub_reord,'col');
        % normalization
        r0_blocks=kernel_col/sum(kernel_col); l0_blocks=fcn_left_kernel(K_sp_sub_reord,r0_blocks,dim_matr);
        % stat sol
        stat_sol_submatr_blocks=r0_blocks*l0_blocks*x0(submatrix_inds);
        stat_sol_blocks(submatrix_inds)=stat_sol_submatr_blocks;
        term_verts_cell{counter}=submatrix_inds;
        
    else
    [vert_topol_sort,term_cycles_ind,~,~,term_cycle_bounds]=fcn_metagraph_scc(A_sparse_sub);
    % [vert_topol_sort,term_cycles_ind,A_metagraph,scc_cell,term_cycle_bounds]=fcn_metagraph_scc(A_orig)
    A_sparse_sub_reordered_terminal=A_sparse_sub(vert_topol_sort,vert_topol_sort);
    K_sp_sub_reord = (A_sparse_sub_reordered_terminal' - speye(dim_matr,dim_matr))*sum(transition_rates_table(:));

    % if cycles are non-terminal, stat sol can be calculated by block inversion, sames as for acyclic graphs
    if isempty(term_cycles_ind)
        % here make sure if 'vert_topol_sort' is the right ordering...
        stat_sol_submatr_blocks=fcn_block_inversion(K_sp_sub_reord,vert_topol_sort,x0,submatrix_inds);
        stat_sol_blocks(submatrix_inds(vert_topol_sort))=stat_sol_submatr_blocks;
        term_verts_cell{counter}=num2cell(submatrix_inds(vert_topol_sort(diag(K_sp_sub_reord)==0)));
    else % if there are terminal cycles, stat sol calc a bit more complicated
         % need to identify terminal cycles, for corresponding columns of
         % kernel we'll need to calculate adjugate matrix
         %
         % probably we don't want it in symbolic form, but just in case
         if isa(K_sp_sub_reord,'double')
            r_null_cycles=zeros(dim_matr,numel(term_cycle_bounds));
         else
            r_null_cycles=sym(zeros(dim_matr,numel(term_cycle_bounds)));
         end
        for k=1:numel(term_cycle_bounds)
            cycle_inds=term_cycle_bounds{k}(1):term_cycle_bounds{k}(end);
            % calc kernel of SCC
            scc_cycle=K_sp_sub_reord(cycle_inds,cycle_inds); % digraph(A_term_cycle(cycle_inds,cycle_inds),'omitselfloops');
            % adjugate matrix -> kernel
            n=numel(cycle_inds); kernel_col=((-1)^(n-1))*fcn_adjug_matrix(scc_cycle,'col');
            % kernel_col_norm=kernel_col/sum(kernel_col);
            r_null_cycles(cycle_inds,k) = kernel_col/sum(kernel_col);
        end
        % if there are single-vertex terminal states too
        if sum(ismember(diag(K_sp_sub_reord),0))>0
            n_terminal=find(ismember(diag(K_sp_sub_reord),0))'; 
            r_null_single_vert = zeros(dim_matr,numel(n_terminal)); 
            r_null_single_vert(sub2ind(size(r_null_single_vert),n_terminal, 1:numel(n_terminal)) )=1;
            % does the order of columns in the kernel matter? I think not, if l0_blocks consistent w r0_blocks
            r0_blocks=[r_null_cycles r_null_single_vert];
        else
            r0_blocks=r_null_cycles;
        end
        % calculate kernel
        l0_blocks=fcn_left_kernel(K_sp_sub_reord,r0_blocks,dim_matr);
        % stat sol
        stat_sol_submatr_blocks=r0_blocks*l0_blocks*x0(submatrix_inds(vert_topol_sort));
        stat_sol_blocks(submatrix_inds(vert_topol_sort))=stat_sol_submatr_blocks;
        % cell of terminal vertices
        [row,col]=find(r0_blocks);
        pre_term_verts_cell=cell(1,numel(unique(col)));
        for k=1:length(pre_term_verts_cell)
            pre_term_verts_cell{k}=submatrix_inds(vert_topol_sort(row(col==k)));
        end
        term_verts_cell{counter}=pre_term_verts_cell;
        
        % term_verts_cell{counter}=submatrix_inds(vert_topol_sort(diag(K_sp_sub_reord)==0));
        
    end % if gate abt terminal cycles

    end % if whole subgraph is one SCC
    
end % end of if gate abt cycles

end % end of if: are there nonzero states in this subgraph?

end % end of for loop going thru disconnected subgraphs

% if there is only one subgraph (no disconnected parts, then cell should have only one layer)
if counter==1
    term_verts_cell=term_verts_cell{1}; cell_subgraphs=cell_subgraphs{1};
end