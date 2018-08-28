% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

classdef Node
    % the Node class contains flow information for each in a network as
    % well as a position for use in plotting the network graphically
    %   Detailed explanation goes here
    
    properties
        inFlow
        outFlow
        outEdges
        inEdges
        X
        Y
        
    end
    
    methods
        function this = Node(f_in, f_out, e_out, e_in, X, Y)
            this.inFlow = f_in;
            this.outFlow = f_out;
            this.outEdges = e_out;
            this.inEdges = e_in;
            this.X = X;
            this.Y = Y;
        end
        
    end
    
end

