classdef Shock
    %SHOCK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M1
        M2
        p1
        p2
        T1
        T2
        rho1
        rho2
        c1
        c2
        theta
        beta
        direction
        type
    end
    
    methods
        function obj = Shock(direction, varargin)
            %SHOCK Class for representing a shock in the shock-fixed frame
            %   Contains all members and methods for representing and
            %   solving a shock in the shock-fixed frame. Note: if in the
            %   lab frame u1=0, then M1=Ms. Generally Ms=u_s/c1,
            %   u_s=|u_shock-u1| (reference frame moving at u_shock). If
            %   theta/beta are unknown, set to -1. The shock type
            %   (weak/strong) only matters when solving for beta.
            if direction == 'R'
                obj.direction = 1;
            else
                obj.direction = -1;
            end
            
            p = inputParser;
            p.CaseSensitive = false;
            addParameter(p,'M1',-1);
            addParameter(p,'M2',-1);
            addParameter(p,'p1',-1);
            addParameter(p,'p2',-1);
            addParameter(p,'T1',-1);
            addParameter(p,'T2',-1);
            addParameter(p,'rho1',-1);
            addParameter(p,'rho2',-1);
            addParameter(p,'c1',-1);
            addParameter(p,'c2',-1);
            addParameter(p,'theta',0);
            addParameter(p,'beta',pi/2);
            addParameter(p,'type','weak');
            parse(p,varargin{:});

            props = properties(obj);
            props = props(not(cellfun(@(x)isequal(x,'direction'),props)));
            
            for iprop = 1:length(props)
                prop = props{iprop};
                obj.(prop) = p.Results.(prop);
            end
        end

        function c = speed_sound(obj, gamma, state)
            if state == 1
                c = sqrt(gamma.*obj.p1./obj.rho1);
            else
                c = sqrt(gamma.*obj.p2./obj.rho2);
            end
        end

        function p = ideal_gas(obj, R, state)
            if state == 1
                p = obj.rho1.*R.*obj.T1;
            else
                p = obj.rho2.*R.*obj.T2;
            end
        end

        function p2p1 = p_ratio(obj, gamma)
            p2p1 = ShockJump.p_ratio(obj.M1,gamma,obj.beta);
        end

        function c2c1 = c_ratio(obj, gamma)
            c2c1 = ShockJump.c_ratio(obj.M1, gamma, obj.beta);
        end

        function h2h1 = h_ratio(obj, gamma)
            h2h1 = ShockJump.h_ratio(obj.M1, gamma, obj.beta);
        end

        function rho2rho1 = rho_ratio(obj, gamma)
            rho2rho1 = ShockJump.rho_ratio(obj.M1, gamma, obj.beta);
        end

        function mach2 = get_M2(obj, gamma)
            mach2 = ShockJump.M2(obj.M1, gamma, obj.beta, obj.theta);
        end

        function t = get_theta(obj, gamma)
            t = ShockJump.theta(obj.M1, gamma, obj.beta);
        end

        function duc1 = du_c1(obj, gamma)
            duc1 = obj.direction.*ShockJump.du_c1(obj.M1, gamma);
        end

        function dpp1 = dp_p1(obj, gamma)
            dpp1 = obj.direction.*ShockJump.dp_p1(obj.du_c1(gamma),gamma);
        end

        function b = get_beta(obj, gamma)
            if strcmp(obj.type, 'weak')
                b0 = asin(1./obj.M1);
            else
                b0 = pi/2;
            end
            fun = @(x)(ShockJump.theta(obj.M1, gamma, x)-obj.theta);
            b = fsolve(fun,b0);
        end

        % function [] = SolveShock(obj, gamma, R)
        %     %SHOCK Solve for unknown shock properties
        %     %   Numerically solves for all unknown shock properties 
        %     %   using Matlab's symbolic toolbox.
        %     %   CURRENTLY NOT WORKING
        % 
        %     if strcmp(obj.type, 'weak')
        %         b0 = asin(1./obj.M1);
        %     else
        %         b0 = pi/2;
        %     end
        % 
        %     props = properties(obj);
        %     unknowns = {};
        %     % Identify all unknown properties
        %     for iprop = 1:length(props)
        %         prop = props{iprop};
        %         if obj.(prop) < 0
        %             obj.(prop) = sym(prop);
        %             assume(obj.(prop), 'real')
        %             unknowns = {unknowns{:} prop};
        %         end
        %     end
        % 
        %     % Sort alphabetically so can match up with output of symvars
        %     unknowns = sort(unknowns);
        % 
        %     % Throw every existing equality at the symbolic solver and hope
        %     % it figures something out
        %     eqns = [
        %         obj.p1 == obj.ideal_gas(R, 1), ...
        %         obj.c1 == obj.speed_sound(gamma, 1), ...
        %         obj.M2 == obj.get_M2(gamma), ...
        %         obj.theta == obj.get_theta(gamma), ...
        %         obj.p2 == obj.p1.*obj.p_ratio(gamma), ...
        %         obj.rho2 == obj.rho1.*obj.rho_ratio(gamma), ...
        %         obj.T2 == obj.T1.*obj.h_ratio(gamma), ...
        %         obj.c2 == obj.c1.*obj.c_ratio(gamma)];
        %     eqns(isAlways(eqns, "Unknown","false")) = [];
        %     x0 = zeros(length(symvar(eqns)),1);
        %     ibeta = find(ismember(unknowns, 'beta'));
        %     if ibeta
        %         x0(ibeta) = b0;
        %     end
        %     soln = vpasolve(eqns, symvar(eqns),x0);
        %     disp(soln);
        % 
        %     % Assign solutions to correct properties
        %     for iprop = 1:length(unknowns)
        %         prop = unknowns{iprop};
        %         obj.(prop) = double(soln.(prop));
        %     end
        % end
        
        function rShock = reflect(obj, gamma, fixed)
            %REFLECT Get reflected shock
            %   Returns a shock object containing the reflected shock,
            %   passing all properties corresponding to its upstream state
            %   (may be used to calculate beta_2 and downstream
            %   properties). The shock type (weak/strong) of the reflected
            %   shock is inherited from the incident shock.
            if fixed
                M = obj.get_M2(gamma);
            else
                syms M_R;
                assume(M_R, 'real');
                M = double(solve((M_R^2-1)/M_R == (1/obj.c_ratio(gamma))*(obj.M1^2-1)/obj.M1, M_R>0));
            end
            rShock = Shock(-obj.direction, ...
                'theta', obj.theta, ...
                'M1', M, ...
                'p1', obj.p1.*obj.p_ratio(gamma), ...
                'T1', obj.T1.*obj.h_ratio(gamma), ...
                'c1', obj.c1.*obj.c_ratio(gamma), ...
                'rho1', obj.rho1.*obj.rho_ratio(gamma), ...
                'type', obj.type);

            % Negative properties are unknown, get rid of them
            props = properties(rShock);
            for iprop = 1:length(props)
                prop = props{iprop};
                if rShock.(prop) < 0
                    rShock.(prop) = -1;
                end
            end
            
            % Enable if I ever get SolveShock working
            % rShock.SolveShock(gamma, R);
        end

        function [] = polar(obj, gamma, reflect)
            %POLAR Plot the shock polar
            %   Plots the shock polar. If reflect is true, also plots the
            %   reflected shock polar.
            hold on;
            
            B1_array = linspace(asin(1/obj.M1), pi/2, 1e3);
            theta1_array = ShockJump.theta(obj.M1,gamma,B1_array);
            p211_array = ShockJump.p_ratio(obj.M1, gamma, B1_array);
            plot(rad2deg([flip(-theta1_array) theta1_array]), ...
                [flip(p211_array) p211_array], 'Color', 'blue');
            
            if reflect
                ref = obj.reflect(gamma, true);
                B2_array = linspace(asin(1/ref.M1), pi/2, 1e3);
                theta2_array = ShockJump.theta(ref.M1, gamma, B2_array);
                p212_array = ShockJump.p_ratio(ref.M1, gamma, B2_array);
                p211 = ShockJump.p_ratio(obj.M1, gamma, obj.beta);
                plot(rad2deg([flip(-theta2_array+ref.theta) theta2_array+ref.theta]), ...
                    [flip(p211.*p212_array) p211.*p212_array], 'Color', 'red');
                xline(0, 'k--')
                legend('Incident Shock', 'Reflected Shock', '', 'Location','northwest');
            end

            xlabel('$\theta$ $[\deg]$', 'Interpreter','latex');
            ylabel('$\frac{p}{p_1}$', 'Interpreter','latex', 'Rotation',0);
            title('Shock Polar')
            
            hold off;
        end
    end
end

