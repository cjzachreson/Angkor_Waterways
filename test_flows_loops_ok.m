% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)


function [ edges_out, nodes_out ] = test_flows_loops_ok(adj_mat, edges, nodes, source_node, sink_node, f_tot )


        in_nodes = source_node; %nodes providing input (sending flow)
        
        while abs(nodes{sink_node}.inFlow - f_tot) / f_tot > (10^(-12))
            
            % loop through nodes sending flow (current_nodes)
            
            out_nodes = [];
            
            for i = in_nodes
                out_nodes = [out_nodes, adj_mat.m{i}];
                % partition flow between output edges
                out_edges = nodes{i}.outEdges;
                
                C_tot = 0;
                
                for eij = out_edges
                    C_tot = C_tot + edges{eij}.capacity;
                end
                
                tot_flow_out = 0;
                
                for eij = out_edges
                    % find intitial flow
                    edges{eij} = ...
                        edges{eij}.find_flow(edges{eij}, nodes{i}.outFlow, C_tot);
                    
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
            
            %abs(nodes{sink_node}.inFlow - f_tot) / f_tot
            
        end

    edges_out = edges;
    nodes_out = nodes;
        
end

