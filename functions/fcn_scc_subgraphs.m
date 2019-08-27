function stg_sorting_cell=fcn_scc_subgraphs(A_sparse,x0)

disp('identifying SCCs')
subnetws=conncomp(digraph(A_sparse,'omitselfloops'),'Type','weak');
num_subnets=length(unique(subnetws)); cell_subgraphs=cell(num_subnets,1); scc_submat_cell=cell(num_subnets,1);

disp('identifying SCCs in subgraphs')
for counter=1:num_subnets
    submatrix_inds=find(subnetws==counter); cell_subgraphs{counter}=submatrix_inds;
    A_sparse_sub=A_sparse(submatrix_inds,submatrix_inds);
    scc_submat_cell{counter}=conncomp(digraph(A_sparse_sub,'omitselfloops'),'Type','strong');
end

nonempty_subgraphs=find(arrayfun(@(x) sum(x0(subnetws==x))>0, 1:numel(scc_submat_cell) ));
sorted_vertices_cell=cell(numel(nonempty_subgraphs),1);
cyclic_sorted_subgraphs_cell=cell(numel(nonempty_subgraphs),1);

counter=0;
disp('Sorting nonempty subgraphs')
for k=nonempty_subgraphs
    counter=counter+1;
    A_sparse_sub=A_sparse(subnetws==k,subnetws==k); 
    % if all SCCs single vertex
    if numel( unique(scc_submat_cell{k}) )==size(A_sparse_sub,1)
        sorted_vertices_cell{counter}=toposort(digraph(A_sparse_sub,'omitselfloops'));
    else % there are multi-vertex SCCs (cycles)
   disp('cycles in STG')
    % if entire graph is one connected component, no reordering needed
    if numel(unique(scc_submat_cell{k}))==1
        sorted_vertices_cell{counter}=find(subnetws==k);
%         K_sp_sub_reord = (A_sparse_sub' - speye(dim_matr,dim_matr))*sum(transition_rates_table(:));
%         kernel_col=((-1)^(dim_matr-1))*fcn_adjug_matrix(K_sp_sub_reord,'col');
%         % normalization
%         r0_blocks=kernel_col'/sum(kernel_col); l0_blocks=fcn_left_kernel(K_sp_sub_reord,r0_blocks,dim_matr);
%         % stat sol
%         stat_sol_submatr_blocks=r0_blocks*l0_blocks*x0(submatrix_inds);
%         stat_sol_blocks(submatrix_inds)=stat_sol_submatr_blocks;
%         term_verts_cell{counter}=submatrix_inds;
    else
        [vert_topol_sort,term_cycles_ind,~,~,term_cycle_bounds]=fcn_metagraph_scc(A_sparse_sub);
        disp(strcat('cycle size',num2str(numel(vert_topol_sort))))
        cyclic_sorted_subgraphs_cell{counter}={vert_topol_sort,term_cycles_ind,term_cycle_bounds};
    end
        
    end
end

stg_sorting_cell = {subnetws,scc_submat_cell,nonempty_subgraphs,sorted_vertices_cell,cyclic_sorted_subgraphs_cell};