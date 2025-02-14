classdef Fan
    %FAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        u1
        u2
        p1
        p2
        T1
        T2
        rho1
        rho2
        c1
        c2
        direction
    end
    
    methods
        function obj = Fan(direction,varargin)
            %FAN Construct an instance of this class
            %   Detailed explanation goes here
            if direction == 'R'
                obj.direction = -1;
            else
                obj.direction = 1;
            end
            p = inputParser;
            addParameter(p,'u1',-1);
            addParameter(p,'u2',-1);
            addParameter(p,'p1',-1);
            addParameter(p,'p2',-1);
            addParameter(p,'T1',-1);
            addParameter(p,'T2',-1);
            addParameter(p,'rho1',-1);
            addParameter(p,'rho2',-1);
            addParameter(p,'c1',-1);
            addParameter(p,'c2',-1);
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
            du_c1 = obj.direction.*(obj.u2-obj.u1)./obj.c1;
            p2p1 = UnsteadyIsentropic.p_ratio(du_c1, gamma);
        end

        function c2c1 = c_ratio(obj, gamma)
            du_c1 = obj.direction.*(obj.u2-obj.u1)./obj.c1;
            c2c1 = UnsteadyIsentropic.c_ratio(du_c1, gamma);
        end

        function h2h1 = h_ratio(obj, gamma)
            du_c1 = obj.direction.*(obj.u2-obj.u1)./obj.c1;
            h2h1 = UnsteadyIsentropic.h_ratio(du_c1, gamma);
        end

        function rho2rho1 = rho_ratio(obj, gamma)
            du_c1 = obj.direction.*(obj.u2-obj.u1)./obj.c1;
            rho2rho1 = UnsteadyIsentropic.rho_ratio(du_c1, gamma);
        end
    end
end

