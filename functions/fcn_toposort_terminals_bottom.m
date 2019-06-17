function [A_sparse_sub_reordered_terminal,sorted_vertices_terminal_bottom]=fcn_toposort_terminals_bottom(A_orig)

if numel(conncomp(digraph(A_orig),'Type','strong'))==size(A_orig,1)

terminal_nodes=find(diag(A_orig)==1)';
% source nodes: source_nodes=find(sum(A_orig - diag(diag(A_orig)) )==0); 
% topological sorting of vertices, order descending with vertices: toposort(G,'Order','stable')
sorted_vertices = toposort(digraph(A_orig,'omitselfloops'));

% this is a consistent ordering but terminals are not necessarily in lower right corner of matrix
A_orig_reordered=A_orig(sorted_vertices,sorted_vertices);
% but we want to have terminal states at the bottom
terminal_indices = find(ismember(sorted_vertices,terminal_nodes)); terminals_rem_inds = find(~ismember(sorted_vertices,terminal_nodes));
sorted_vertices_terminal_bottom = [sorted_vertices(~ismember(sorted_vertices, terminal_nodes)) sorted_vertices(terminal_indices)];
% vert_topol_sort_term_bottom=[terminals_rem_inds terminal_indices];

A_sparse_sub_reordered_terminal = A_orig_reordered([terminals_rem_inds terminal_indices],[terminals_rem_inds terminal_indices]);

else
   disp('contains cycles!!') 
end