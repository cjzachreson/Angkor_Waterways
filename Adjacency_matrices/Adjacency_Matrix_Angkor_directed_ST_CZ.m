% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Adjacency_Matrix_Angkor_directed_ST_CZ
    % adjacency matrix for a directed network, each row index corresponds
    % to a node ID while each column entry is a node to which flow is
    % directed from the node corresponding to the row index
   
    
    properties
        m %sparse adjacency matrix specifying connected nodes
        v %values of edges between nodes specified by m
        N
        name
    end
    
    methods
        
        function this = Adjacency_Matrix_Angkor_directed_ST_CZ(root)
           
            % initializes the Angkor network by reading it in from an
            % external source and translating the full adjacency matrix to a
            % sparse one. 
            
            this.name = 'ankgor_CZ_dir_ST';
            
            full_matrix = dlmread(root);

            this.m = {};
            this.v = {};
            
            % first, identify source and sink nodes
            % for now, I assume there are more sinks than sources
            % currently I assume the matrix dimension denotes the edge
            % direction ie row->column or column->row
            
            size_m = size(full_matrix, 1);
            
            %assuming matrix does not currently include super-source and
            %super-sink. 
            this.N = size_m;
            
%            topology = sign(full_matrix);g = biograph(topology);view(g);

            %now go through all the rows (output nodes) assigning connections and
            %values. 
            
            for i = 1:size_m-1
                this.m{i} = find(full_matrix(i, :));
                this.v{i} = full_matrix(i, this.m{i});
            end
            
            this.m{this.N} = [];
            this.v{this.N} = [];
            
        end
        
    end
    
end














