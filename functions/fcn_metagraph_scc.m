function [vert_topol_sort,term_cycles_ind,A_metagraph,scc_cell,term_cycle_bounds]=fcn_metagraph_scc(A_sparse_sub)

matr_size=size(A_sparse_sub); scc_list=conncomp(digraph(A_sparse_sub,'omitselfloops'),'Type','strong');
[num_verts_per_scc,scc_memb_per_vert]=histc(scc_list,unique(scc_list));
scc_cell=conncomp(digraph(A_sparse_sub),'OutputForm','cell');

numel_scc=numel(num_verts_per_scc);
A_metagraph = zeros(numel_scc,numel_scc); size_A_metagr=size(A_metagraph);
% we need to convert indices
[row,col]=ind2sub(matr_size,find((A_sparse_sub>0)- diag(diag(A_sparse_sub)) ) );
% remove those that are within SCCs (since we remove diag elements, these belong to the same SCCs)
rel_edges=find(scc_memb_per_vert(row)~=scc_memb_per_vert(col));
metagraph_edges=scc_memb_per_vert([row(rel_edges),col(rel_edges)]); 
A_metagraph(sub2ind(size_A_metagr,metagraph_edges(:,1),metagraph_edges(:,2))) = A_sparse_sub(sub2ind(matr_size,row(rel_edges),col(rel_edges)));
% are terminal vertices in the lower right block of matrix?

metagraph_ordering=toposort(digraph(A_metagraph));
terminal_scc_ind=find(sum(A_metagraph,2)==0)'; terminal_scc_pos=ismember(metagraph_ordering,terminal_scc_ind);
nonterm_scc_num=numel_scc-numel(terminal_scc_ind);
% identify terminal cycles
term_cycles_ind=intersect(find(cellfun(@(x) numel(x),scc_cell)>1),terminal_scc_ind);

if sum(~(terminal_scc_pos>nonterm_scc_num))>0
    nonterm_scc_inds = ~ismember(metagraph_ordering,terminal_scc_ind);
    metagraph_ordering_terminal_bottom = [metagraph_ordering(nonterm_scc_inds) metagraph_ordering(terminal_scc_pos)];
end

if ~isempty(term_cycles_ind)
    scc_cell_reordered=scc_cell(metagraph_ordering_terminal_bottom);
    % index of cells containing term cycles after reordering
    term_cycles_ind=find(ismember(metagraph_ordering_terminal_bottom,term_cycles_ind));
    % we need a cell of the indices of vertices within these
    scc_cell_reordered_cumsum=cumsum(cellfun(@(x) numel(x),scc_cell_reordered));
    scc_cell_reordered_lengths=cellfun(@(x) numel(x),scc_cell_reordered);
    cycle_first_verts=scc_cell_reordered_cumsum(term_cycles_ind) - scc_cell_reordered_lengths(term_cycles_ind)+1;
    cycle_last_verts=scc_cell_reordered_cumsum(term_cycles_ind);
    term_cycle_bounds = num2cell([cycle_first_verts; cycle_last_verts]',2);
    
else
    term_cycles_ind=[];
    term_cycle_bounds=[];
end

% reordered original vertices
vert_topol_sort=cell2mat(scc_cell(metagraph_ordering_terminal_bottom));
