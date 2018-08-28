% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

function [ edges, nodes ] = initialize_network_valued( adj_mat, alpha, beta )

        edges = {};
        nodes = {};
        edge_ID = 0;
        
        for i = 1:adj_mat.N;
            
            out_nodes = adj_mat.m{i};
            edge_values = adj_mat.v{i};
            
            if numel(out_nodes) ~= numel(edge_values)
                error(['unvalued edge detected in valued network'])
            end
            
            n_out = numel(out_nodes);
            
            n_prev = edge_ID;
            
            n_current = n_prev + numel(out_nodes);
            
            out_edges = n_prev + 1 : 1 : n_current;
            % may set X and Y manually
            nodes{i} = Node(0, 0, out_edges, [], [], []);
            
            for j = 1:n_out
                
                edge_ID = edge_ID + 1;
                
                competitors = out_edges(out_edges ~= edge_ID);
                
                edges{edge_ID} = Edge(edge_ID, i, out_nodes(j), edge_values(j), alpha, beta, competitors);
                
            end
            
        end
        
        for i = 1:numel(edges)
            nodes{edges{i}.outNode}.inEdges = [nodes{edges{i}.outNode}.inEdges, i];
        end
        
   
end

