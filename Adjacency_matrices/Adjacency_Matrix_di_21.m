% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Adjacency_Matrix_di_21
    % adjacency matrix for a directed network, each row index corresponds
    % to a node ID while each column entry is a node to which flow is
    % directed from the node corresponding to the row index
   
    
    properties
        m
        N
    end
    
    methods
        
        function this = Adjacency_Matrix_di_21()
           
            %initializes a 9-member pyramidal network with supersource and
            %supersink
            this.N = 21;
            
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
            this.m{10} = [14, 15];
            this.m{11} = [15, 16];
            this.m{12} = [16, 17];
            this.m{13} = [17];
            
            this.m{14} = [18];
            this.m{15} = [18, 19];
            this.m{16} = [19, 20];
            this.m{17} = [20];
            
            this.m{18} = [21];
            this.m{19} = [21];
            this.m{20} = [21];
            
            
            this.m{21} = [];
            
        end
        
    end
    
end

