function x0=fcn_define_initial_states(initial_on_nodes,dom_prob,nodes,distrib_type,plot_flag)

n_nodes=numel(nodes); truth_table_inputs=rem(floor([0:((2^n_nodes)-1)].'*pow2(0:-1:-n_nodes+1)),2);
% define initial values
x0=zeros(1,2^n_nodes)'; 
% defining a dominant initial state (eg. dom_prob=0.8, ie. 80% probability)
initial_on_nodes_inds=find(sum(cell2mat(arrayfun(@(x) strcmp(nodes,initial_on_nodes{x} ), 1:numel(initial_on_nodes),'un',0)')));
% dom_prob=0.8; 
initial_state=zeros(1,numel(nodes)); initial_state(initial_on_nodes_inds)=1;
x0(ismember(truth_table_inputs,initial_state,'rows'))=dom_prob; 

% take those states where selected variables have a value of 1, and we want them to have a nonzero probability
if strcmp(distrib_type,'restrict')
sel_states=all(truth_table_inputs(:,initial_on_nodes_inds)');
% create a vector of random probabilities for these states, with a sum of (1-dom_prob)
rand_vect=abs(rand(sum(sel_states)-1,1)); rand_vect=rand_vect/( sum(rand_vect)/(1-dom_prob) );
x0(~ismember(truth_table_inputs,initial_state,'rows') & sel_states') = rand_vect;
elseif strcmp(distrib_type,'broad')
    sel_states = ~ismember(truth_table_inputs,initial_state,'rows')';
    % create a vector of random probabilities for these states, with a sum of (1-dom_prob)
    rand_vect=abs(rand(sum(sel_states),1)); rand_vect=rand_vect/( sum(rand_vect)/(1-dom_prob) );
    x0(sel_states) = rand_vect;
else
   error('<distrib_type> should be ''restrict'' or ''broad''.') 
end

% rounding precision
n_prec=3;
if round(sum(x0),n_prec)==1
    disp('sum(x0)=1, OK.')
else
     disp('sum(x0)~=1, something wrong!')
end

if ~isempty(plot_flag)
bar(x0); set(gca,'yscale','log'); xlim([1 2^n_nodes]); ylim([(1-dom_prob)/2^n_nodes 1])
% subplot(2,1,2); x0=fcn_define_initial_states(initial_on_nodes,dom_prob,nodes,'broad'); bar(x0); xlim([1 2^13]);set(gca,'yscale','log'); ylim([(1-dom_prob)/2^n_nodes 1])
end