% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Adjacency_Matrix_simple_4_val
    % adjacency matrix for a directed network, each row index corresponds
    % to a node ID while each column entry is a node to which flow is
    % directed from the node corresponding to the row index
   
    
    properties
        m
        v
        N
        name
    end
    
    methods
        
        function this = Adjacency_Matrix_simple_4_val()
           this.name = 'simple_4_val';
            %initializes a 9-member pyramidal network with supersource and
            %supersink
            this.N = 4;
            this.m = {}; this.v = {};
            this.m{1} = [2, 3]; this.v{1} = [1, 1];
            this.m{2} = [4]; this.v{2} = [1];
            this.m{3} = [4]; this.v{3} = [1];
            this.m{4} = []; this.v{4} = [];

            
        end
        
    end
    
end

