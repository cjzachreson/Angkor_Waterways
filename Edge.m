
% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)



classdef Edge
    % Edge class contians information about the flow through an edge, the
    % incoming node and the outgoing node, the other edges attached to the
    % same incoming node, the errodabilty of the edge (for irrigation
    % network simulation), and the efficiency of the edge
    
    % methods include calcuation of flow based on the incoming flow and the
    % efficiencies of competing edges, and the dynamic function that
    % describes how efficiency changes as a function of flow and efficiency
    
    %   Detailed explanation goes here
    
    properties
        ID
        inNode
        outNode
        initialization_weight
        capacity
        initial_capacity;
        flow
        initial_flow
        tot_flow_t
        alpha
        beta
        competitors
        erosion_threshold
        sedimentation_threshold
        damage
        velocity
        initial_velocity
        
    end
    
    methods (Static)
        
        function this = Edge(id, inNode, outNode, w, alpha, beta, competitors)
            this.ID = id;
            this.inNode = inNode;
            this.outNode = outNode;
            this.initialization_weight = w;
            this.alpha = alpha;
            this.beta = beta;
            this.competitors = competitors;
            this.tot_flow_t = 0;
        end
        
        function this = initialize_flow(this, f_in, w_tot)
            this.flow = f_in *(this.initialization_weight / w_tot);
        end
        
        function this = find_flow(this, f_in, C_tot)
            this.flow = f_in * (this.capacity / C_tot);
        end
        
        function this = set_thresholds(this)
            this.erosion_threshold = this.velocity / this.alpha; %velocity
            this.sedimentation_threshold = this.beta * this.velocity;
        end
        
        function this = set_velocity(this)
            this.velocity = this.tot_flow_t / this.capacity;
        end
        
        function this = update_capacity(this) %veocity = flow / efficiency
            
            if this.velocity > this.erosion_threshold
                this.capacity = ...
                    this.capacity + this.alpha * this.flow * ...
                    (this.velocity / this.erosion_threshold - 1);
            end
            
            if this.velocity < this.sedimentation_threshold
                this.capacity = ...
                    this.capacity - this.beta * this.flow * ...
                    (this.sedimentation_threshold / this.velocity - 1);
                
                % need to prevent numerical errors when efficiency is near zero                
                small_num = 1 * 10^(-12);
                
                if this.capacity < 0
                    error(['negative capacity on edge ' num2str(this.ID)])
                end
                
                if this.capacity < small_num
                    this.capacity = small_num;
                end
                
            end
            
        end
        
    end
    
end

