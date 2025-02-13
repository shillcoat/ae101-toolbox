classdef RiemannProblem
    %RIEMANNPROBLEM Plotting and solution of the canonical Riemann problem.
    %   Detailed explanation goes here
    
    properties
        p_L
        u_L
        T_L
        c_L
        g
    end
    
    methods
        function obj = RiemannProblem(pL, uL, TL, cL, gamma)
            %RIEMANNPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            obj.p_L = pL;
            obj.u_L = uL;
            obj.T_L = TL;
            obj.c_L = cL;
            obj.g = gamma;
        end
        
        function p_pL = S_neg(obj,u_cL)
            %S_NEG Summary of this method goes here
            %   Detailed explanation goes here
            p_pL = 1 + ShockJump.dp_p1(obj.u_L./obj.c_L-u_cL, obj.g);
        end

        function p_pL = S_hat_pos(obj,u_cL)
            %S_HAT_POS Summary of this method goes here
            %   Detailed explanation goes here
            p_pL = 1 - obj.g.*(obj.u_L/obj.c_L-u_cL).*...
                (-(obj.g+1)/4.*(obj.u_L./obj.c_L-u_cL)+...
                sqrt(1+(obj.g+1).^2./16*(obj.u_L./obj.c_L-u_cL).^2));
        end

        function p_pL = R_neg(obj,u_cL)
            %R_NEG Summary of this method goes here
            %   Detailed explanation goes here
            p_pL = UnsteadyIsentropic.p_ratio(obj.u_L./obj.c_L-u_cL,obj.g);
        end

        function p_pL = R_hat_pos(obj,u_cL)
            %R_HAT_POS Summary of this method goes here
            %   Detailed explanation goes here
            p_pL = UnsteadyIsentropic.p_ratio(u_cL-obj.u_L./obj.c_L,obj.g);
        end

        function [] = plot(obj, xrange)
            % PLOT Summary of this method goes here
            %   Detailed explanation goes here
            x = linspace(xrange(1), xrange(2), 1e3);
            L_ind = find(x>=obj.u_L/obj.c_L);
            L_ind = L_ind(:,1);

            plot(x(1:L_ind),obj.S_neg(x(1:L_ind)),...
                'LineStyle', '-','Color','b','LineWidth',2);
            hold on;
            plot(x(L_ind:end),obj.R_neg(x(L_ind:end)),...
                'LineStyle', '-','Color','k','LineWidth',2);
            plot(x(1:L_ind),obj.S_hat_pos(x(1:L_ind)),...
                'LineStyle','--','Color','b','LineWidth',2);
            plot(x(L_ind:end),obj.R_hat_pos(x(L_ind:end)),...
                'LineStyle','--','Color','k','LineWidth',2);
            legend(["$S^-$", "$R^-$", "$\hat{S}^+$", "$\hat{R}^+$"], ...
                "Interpreter", "latex", "Location", "north", ...
                'FontSize', 12);
            xlabel('$u/c_L$', 'Interpreter','latex','FontSize',14);
            ylabel('$p/p_L$', 'Interpreter','latex','FontSize',14);
            hold off;
        end
    end
end

