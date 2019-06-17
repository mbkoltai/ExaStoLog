function stat_sol_submatr_blocks=fcn_block_inversion(K_sp_sub_reord,sorted_vertices_terminal_bottom,x0,submatrix_inds)

% this function calculates kernels and stationary solution if all terminal
% SCCs are single vertices (states)

% construct kernels from matrix blocks
dim_kernel=sum(diag(K_sp_sub_reord)==0); dim_matr=size(K_sp_sub_reord,1); % bc this is no cycles branch of IF gate
colnum_r_null_array=1:dim_kernel; 
term_block_inds = (dim_matr-dim_kernel+1):dim_matr; nonterm_block_inds = 1:dim_matr-dim_kernel;
term_block = speye(dim_kernel,dim_kernel);
% right kernel
r0_blocks=zeros(dim_matr,dim_kernel); r0_blocks(term_block_inds, colnum_r_null_array)=term_block; r0_blocks=sparse(r0_blocks);
% left kernel
l0_blocks = transpose(zeros(size(r0_blocks))); l0_blocks(transpose(~ismember(r0_blocks,0)))=1; 
if isa(r0_blocks,'sym'); l0_blocks=sym(l0_blocks); end
X_block = -r0_blocks(term_block_inds,colnum_r_null_array)*K_sp_sub_reord(term_block_inds,nonterm_block_inds)/(K_sp_sub_reord(nonterm_block_inds,nonterm_block_inds));
l0_blocks(colnum_r_null_array,nonterm_block_inds)=X_block;
% due to reordering of states init cond state inds need to be reordered too
stat_sol_submatr_blocks = r0_blocks*l0_blocks*x0(submatrix_inds(sorted_vertices_terminal_bottom));
