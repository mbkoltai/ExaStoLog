function adj_matrix=fcn_adjug_matrix(A,col_arg)


size_array=1:size(A,1); size_vect=numel(size_array);
if isa(A,'double')
    adj_matrix = zeros(size(A));
elseif isa(A,'sym')
    adj_matrix = sym(zeros(size(A)));
end

if size(A,1)==size(A,2)
    
if isempty(col_arg)
        
for k=1:size_vect
for l = 1:size_vect
    adj_matrix(k,l) = ((-1)^(k+l))*det(A(setdiff(size_array,l),setdiff(size_array,k)));
end
end

if isa(adj_matrix,'sym')
    adj_matrix=simplify(adj_matrix);
end

else % if we don't want all cols, jmust 1st (since they are identical)
    
adj_matrix=adj_matrix(:,1);
l=1;

for k=1:size_vect
    adj_matrix(k) = ((-1)^(k+l))*det(A(setdiff(size_array,l),setdiff(size_array,k)));
end 

end

else % if non-square
    adj_matrix=[];
    disp('non-square matrix!!')
end