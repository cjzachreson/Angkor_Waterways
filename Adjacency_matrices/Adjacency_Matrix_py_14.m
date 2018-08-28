% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Adjacency_Matrix_py_14
    % adjacency matrix for a directed network, each row index corresponds
    % to a node ID while each column entry is a node to which flow is
    % directed from the node corresponding to the row index
   
    
    properties
        m
        N
    end
    
    methods
        
        function this = Adjacency_Matrix_py_14()
           
            %initializes a 9-member pyramidal network with supersource and
            %supersink
            this.N = 14;
            this.m = {};
            
            this.m{1} = [2, 3, 4];
            
            this.m{2} = [5, 6];
            this.m{3} = [6, 7];
            this.m{4} = [7, 8];
            
            this.m{5} = [9, 10];
            this.m{6} = [10, 11];
            this.m{7} = [11, 12];
            this.m{8} = [12, 13];
            
            this.m{9} = [14];
            this.m{10} = [14];
            this.m{11} = [14];
            this.m{12} = [14];
            this.m{13} = [14];
            
            this.m{14} = [];
            
        end
        
    end
    
end

