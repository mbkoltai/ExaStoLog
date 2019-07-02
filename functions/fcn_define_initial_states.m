function x0=fcn_define_initial_states(initial_fixed_nodes,initial_fixed_nodes_vals,dom_prob,nodes,distrib_type,plot_flag)
                                   
n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
% define initial values
x0=zeros(2^n_nodes,1); 
% defining a dominant initial state (eg. dom_prob=0.8, ie. 80% probability)
initial_on_nodes_inds=any(cell2mat(arrayfun(@(x) strcmp(nodes,initial_fixed_nodes{x} ), 1:numel(initial_fixed_nodes),'un',0)'));
% dom_prob=0.8; 
% initial_state=zeros(1,numel(nodes)); initial_state(initial_on_nodes_inds)=1;
inds_condition=ismember(truth_table_inputs(:,initial_on_nodes_inds),initial_fixed_nodes_vals,'rows');
if strfind(distrib_type,'unif')
    x0(inds_condition)=repmat(dom_prob/sum(inds_condition),sum(inds_condition),1); 
    x0(~inds_condition)=repmat((1-dom_prob)/(numel(x0)-sum(inds_condition)),numel(x0)-sum(inds_condition),1); 
elseif strfind(distrib_type,'rand')
    x0(inds_condition)=rand(sum(inds_condition),1); x0=dom_prob*x0/sum(x0);
    x0(~inds_condition)=rand(numel(x0)-sum(inds_condition),1); x0(~inds_condition)=(1-dom_prob)*x0(~inds_condition)/sum(x0(~inds_condition));
else
    error('distrib type shoul be ''uniform'' or ''random''')
end
    
% take those states where selected variables have a value of 1, and we want them to have a nonzero probability
% if strcmp(alloc_type,'restrict')
% sel_states=all(truth_table_inputs(:,initial_on_nodes_inds)');
% % create a vector of random probabilities for these states, with a sum of (1-dom_prob)
%     if strcmp(distrib_type,'random')
%     rand_vect=abs(rand(sum(sel_states)-1,1)); rand_vect=rand_vect/( sum(rand_vect)/(1-dom_prob) );
%     x0(~ismember(truth_table_inputs,initial_state,'rows') & sel_states') = rand_vect; % 
% % other selected states have uniform probability:
%     elseif strfind(distrib_type,'unif')
%     x0(~ismember(truth_table_inputs,initial_state,'rows') & sel_states') = (1-dom_prob)/(sum(sel_states)-1); % 
%     end
% 
% elseif strcmp(alloc_type,'broad')
%     sel_states = ~ismember(truth_table_inputs,initial_state,'rows')';
%     % create a vector of random probabilities for these states, with a sum of (1-dom_prob)
% 	if strcmp(distrib_type,'random')
%         rand_vect=abs(rand(sum(sel_states),1)); rand_vect=rand_vect/( sum(rand_vect)/(1-dom_prob) );
%         x0(sel_states)=rand_vect;
%     elseif strfind(distrib_type,'unif')
%         x0(sel_states) = (1-dom_prob)/(sum(sel_states)-1);
% 	end
%     
% else
%    error('<alloc_types> should be ''restrict'' or ''broad''.') 
% end

% rounding precision
n_prec=3;
if round(sum(x0),n_prec)==1
    disp('sum(x0)=1, OK.')
else
     disp('sum(x0)~=1, something wrong!')
end

if ~isempty(plot_flag)
bar(x0); set(gca,'yscale','log'); xlim([1 2^n_nodes]); % ylim([(1-dom_prob)/2^n_nodes 1])
% subplot(2,1,2); x0=fcn_define_initial_states(initial_on_nodes,dom_prob,nodes,'broad'); bar(x0); xlim([1 2^13]);set(gca,'yscale','log'); ylim([(1-dom_prob)/2^n_nodes 1])
end