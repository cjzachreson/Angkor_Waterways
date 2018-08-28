% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Adjacency_Matrix_py_N_val
    % adjacency matrix for a directed network, each row index corresponds
    % to a node ID while each column entry is a node to which flow is
    % directed from the node corresponding to the row index
   
    
    properties
        m
        N
        v
        name
    end
    
    methods
        
        function this = Adjacency_Matrix_py_N_val(n_levels)
           
            %initializes a n-level pyramidal network with supersource and
            %supersink
            this.N = max(cumsum(1:n_levels)) + 1;
            this.m = {};
            this.name = ['py_' num2str(n_levels) '_levels'];
            
            for i = 1:n_levels
                
                level_nodes = (max(cumsum(0:i-1)) + 1) :  max(cumsum(0:i));
                next_level_nodes = (max(cumsum(0:i)) + 1) :  max(cumsum(0:i+1));
                
                if i ~= n_levels
                
                    for j = 1:numel(level_nodes)
                        this.m{level_nodes(j)} = [next_level_nodes(j), next_level_nodes(j+1)];
                        this.v{level_nodes(j)} = [1, 1];
                    end
                    
                else
                    
                    for j = 1:numel(level_nodes)
                        this.m{level_nodes(j)} = [this.N];
                        this.v{level_nodes(j)} = [1];
                    end
                    
                end
                
            end
            
            this.m{this.N} = [];
            this.v{this.N} = [];
                
                
            
            
%             this.m{1} = [2, 3];

%             this.m{2} = [4, 5];
%             this.m{3} = [5, 6];
%          
%             this.m{4} = [7, 8];
%             this.m{5} = [8, 9];
%             this.m{6} = [9, 10];

%             this.m{7} = [11];
%             this.m{8} = [11];
%             this.m{9} = [11];                
%             this.m{10} = [11];

%             this.m{11} = [];
            
        end
        
    end
    
end

