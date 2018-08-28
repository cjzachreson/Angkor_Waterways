% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Adjacency_Matrix_py_11_val
    % adjacency matrix for a directed network, each row index corresponds
    % to a node ID while each column entry is a node to which flow is
    % directed from the node corresponding to the row index
   
    
    properties
        m
        v
        N
    end
    
    methods
        
        function this = Adjacency_Matrix_py_11_val()
           
            %initializes a 9-member pyramidal network with supersource and
            %supersink
            this.N = 11;
            this.m = {};
            this.m{1} = [2, 3]; this.v{1} = [1, 1];
            this.m{2} = [4, 5]; this.v{2} = [1, 1];
            this.m{3} = [5, 6]; this.v{3} = [1, 1];
            this.m{4} = [7, 8]; this.v{4} = [1, 1];
            this.m{5} = [8, 9]; this.v{5} = [1, 1];
            this.m{6} = [9, 10]; this.v{6} = [1, 1];
            this.m{7} = [11]; this.v{7} = [1];
            this.m{8} = [11]; this.v{8} = [1];
            this.m{9} = [11]; this.v{9} = [1];
            this.m{10} = [11]; this.v{10} = [1];
            this.m{11} = []; this.v{11} = [];
            
        end
        
    end
    
end

