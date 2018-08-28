% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

function [ edges_out, nodes_out ] = initialize_flows_loops_ok(adj_mat, edges, nodes, source_node, sink_node, f_tot )
        % loops converge exponentially towards steady-state flows, this
        % code implements a small-number cuttoff at which point convergance
        % of the flow through the network is established. 

        in_nodes = source_node; %nodes providing input (sending flow)
        
        while abs(nodes{sink_node}.inFlow - f_tot) / f_tot > (10^(-12))
            
            % loop through nodes sending flow (current_nodes)
            
            out_nodes = [];
            
            for i = in_nodes
                out_nodes = [out_nodes, adj_mat.m{i}];
                % partition flow between output edges
                out_edges = nodes{i}.outEdges;
                
                w_tot = 0;
                
                for eij = out_edges
                    w_tot = w_tot + edges{eij}.initialization_weight;
                end
                
                tot_flow_out = 0;
                
                for eij = out_edges
                    % find intitial flow
                    edges{eij} = ...
                        edges{eij}.initialize_flow(edges{eij}, nodes{i}.outFlow, w_tot);
                    
                    tot_flow_out =  tot_flow_out + edges{eij}.flow;
                    edges{eij}.tot_flow_t = edges{eij}.tot_flow_t + edges{eij}.flow;
                    
                    nodes{edges{eij}.outNode}.inFlow =...
                        nodes{edges{eij}.outNode}.inFlow + edges{eij}.flow;
                    
                    nodes{edges{eij}.outNode}.outFlow =...
                        nodes{edges{eij}.outNode}.outFlow + edges{eij}.flow;
                    
                end
                
                nodes{i}.outFlow = nodes{i}.outFlow - tot_flow_out;
                
            end
            
            in_nodes = unique(out_nodes);
            
        end

    edges_out = edges;
    nodes_out = nodes;
        
end

