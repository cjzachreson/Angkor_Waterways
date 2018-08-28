% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Adjacency_Matrix_py_11
    % adjacency matrix for a directed network, each row index corresponds
    % to a node ID while each column entry is a node to which flow is
    % directed from the node corresponding to the row index
   
    
    properties
        m
        N
    end
    
    methods
        
        function this = Adjacency_Matrix_py_11()
           
            %initializes a 9-member pyramidal network with supersource and
            %supersink
            this.N = 11;
            this.m = {};
            this.m{1} = [2, 3];
            this.m{2} = [4, 5];
            this.m{3} = [5, 6];
            this.m{4} = [7, 8];
            this.m{5} = [8, 9];
            this.m{6} = [9, 10];
            this.m{7} = [11];
            this.m{8} = [11];
            this.m{9} = [11];
            this.m{10} = [11];
            this.m{11} = [];
            
        end
        
    end
    
end

