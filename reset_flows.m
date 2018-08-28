% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

function [edges_out, nodes_out] = reset_flows( edges_in, nodes_in )

edges_out = edges_in;
nodes_out = nodes_in;

% reset flows
for n_i = 1:numel(nodes_out)
    nodes_out{n_i}.inFlow = 0;
    nodes_out{n_i}.outFlow = 0;
end

for eij = 1:numel(edges_out)
    edges_out{eij}.tot_flow_t = 0;
    edges_out{eij}.flow = 0;
end


end

