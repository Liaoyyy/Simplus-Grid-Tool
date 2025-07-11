 % This class defines the model of a photovoltaic system.

% Author(s): Wenjie Ning
%


%% Notes
%
% The model is in load convention.
% 
% The model is in admittance form.
%
% dw means the derivative of w

%% Class

classdef Photovoltaic < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        v_od_r;
        v_oq_r;
        P0;
        Q0;
    end
    
    methods
        % constructor
        function obj = Photovoltaic(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            State  = {'v_i','i_i','i_l','v_pv','v_di','v_qi','i_di','i_qi','i_d','i_q','v_d','v_q','i_gd','i_gq','v_o','theta'};  
            Input  = {'v_d','v_q','P0'};
            Output = {'i_d','i_q','w','theta'};
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            % Get parameter
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = 0;
            W0      = obj.Para(11);
            S_air   = obj.Para(12);
            T_air   = obj.Para(13);
            v_pv    = obj.Para(15);
            v_o     = obj.Para(16);
            % Solar para
            T_ref = 25;
            S_ref = 1000;
            I_sc  = 14.880;
            I_m   = 13.88;
            U_m   = 576;
            U_oc  = 708;
            a    = 0.0025;
            b    = 0.00288;
            Rs   = 0.5;
            k    = 0.03;
            C2   = (U_m/U_oc-1)/(log(1-I_m/I_sc));
            C1   = (1-I_m/I_sc)*exp(-U_m/(C2*U_oc));
            T    = T_air + k*S_air;
            dT   = T - T_ref;
            dI   = a * S_air/S_ref * dT+I_sc*(1-S_ref/S_air);
            dU   = -b * dT - Rs * dI;

            % Calculate parameters
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;
            Rf_dc =0.01;
            
            v_gd = V;
            v_gq = 0;
            i_gd = P/V;
            i_gq = -Q/V;
            
            v_gdq = v_gd + 1i*v_gq;
            i_gdq = i_gd + 1i*i_gq;
            v_dq = v_gdq - i_gdq*(Rc + 1i*w*Lc);
            i_cdq = v_dq*(1i*w*Cf);
            i_dq = i_gdq - i_cdq;
            e_dq  = v_dq - i_dq*(Rf + 1i*w*Lf);
            
            i_d = real(i_dq);
            i_q = imag(i_dq);
            i_d_i = -real(e_dq);         
            i_q_i = -imag(e_dq);
            v_d = real(v_dq);
            v_q = imag(v_dq);
            v_d_i = -i_d;              
            v_q_i = -i_q;
            i_gd = real(i_gdq);
            i_gq = imag(i_gdq);
            theta = xi;

            % dc
            i_pv = (( I_sc * (1-C1 * (exp((v_pv*800+dU)/C2/U_oc)-1)))+dI)/40;
            i_l    = i_pv;
            v_i    = i_l; 
            ed_dc = v_pv-i_l*Rf_dc;
            i_i    = ed_dc;

            obj.P0 = P*(-1);
            obj.Q0 = Q*(-1);
            
            v_odq_r = v_dq + (Rov + 1i*Xov)*i_gdq*(-1);
            v_od_r = real(v_odq_r);
            v_oq_r = imag(v_odq_r);
            obj.v_od_r = v_od_r;
            obj.v_oq_r = v_oq_r;
            
            % Get equilibrium
            x_e = [v_i;i_i;i_l;v_pv;v_d_i;v_q_i;i_d_i;i_q_i;i_d;i_q;v_d;v_q;i_gd;i_gq;v_o;theta];
            u_e = [v_gd; v_gq; P];
        end
      
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Set the invertor control type
            Control_type = 1;  %% 1--Closed loop 0--Open loop

            % Get the power PowerFlow values
            V	= obj.PowerFlow(3);
            
            % Get input
            vgd   = u(1);
            vgq   = u(2);

            % Get state
            v_i    = x(1);
            i_i    = x(2);
            i_l    = x(3);
            v_pv   = x(4);
            v_d_i   = x(5);
            v_q_i   = x(6); 
            i_d_i   = x(7);
            i_q_i   = x(8);
            i_d    = x(9);
            i_q    = x(10);
            v_d    = x(11);
            v_q    = x(12);
            i_gd   = x(13);
            i_gq   = x(14);
            v_o    = x(15);
            theta = x(16);
            
            % Get parameters
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = 0;
            N1      = obj.Para(7);        
            R1      = obj.Para(8);            
            xfvdq   = obj.Para(9);
            xfidq   = obj.Para(10);
            W0      = obj.Para(11);
            w0      = W0 ; 
            S_air   = obj.Para(12);
            T_air   = obj.Para(13);
            Co      = obj.Para(14); 
            vr      = obj.Para(15);
            V_dc    = obj.Para(16);
            xfvdc   = obj.Para(17);
            xfidc   = obj.Para(18);

            % Solar para
            T_ref = 25;
            S_ref = 1000;
            I_sc  = 14.880;
            I_m   = 13.88;
            U_m   = 576;
            U_oc  = 708;
            a    = 0.0025;
            b    = 0.00288;
            Rs   = 0.5;
            k    = 0.03;
            C2   = (U_m/U_oc-1)/(log(1-I_m/I_sc));
            C1   = (1-I_m/I_sc)*exp(-U_m/(C2*U_oc));
            T    = T_air + k*S_air;
            dT   = T - T_ref;
            dI   = a * S_air/S_ref * dT+I_sc*(1-S_ref/S_air);
            dU   = -b * dT - Rs * dI;

            % Update paramters
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;

            Rf_dc = 0.01;
            Lf_dc = 0.05/W0;
            Cf_dc = 0.02/W0;

            wv_dc = xfvdc*2*pi; 
            kpv_dc = -Cf_dc*wv_dc* 20;
            kiv_dc = -Cf_dc*wv_dc^2/4* 20;

            wi_dc  = xfidc*2*pi;
            kpi_dc = -Lf_dc*wi_dc;
            kii_dc = -Lf_dc*(wi_dc^2)/4;

            wv1  = xfvdq*2*pi;    % Voltage Controller
            kpv1 = Cf*wv1;
            kiv1 = Cf*wv1^2/4* 20;

            wi1  = xfidq*2*pi;  % Current Controller
            kpi1 = Lf*wi1;
            kii1 = Lf*(wi1^2)/4;
       
  
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                i_pv = (( I_sc * (1-C1 * (exp((v_pv*800+dU)/C2/U_oc)-1)))+dI)/40;
                error_v_pv = vr-v_pv;
                i_r        = kpv_dc*error_v_pv + v_i;
                dv_i       = kiv_dc *error_v_pv;

                error_i_l  = i_r - i_l;
                ed_dc        = kpi_dc*error_i_l + i_i;
                di_i       = kii_dc*error_i_l;

                di_l       = (v_pv-ed_dc-Rf_dc*i_l)/Lf_dc;
                dv_pv      = (i_pv-i_l)/Cf_dc;
                
                v_dr = obj.v_od_r;
                v_qr = obj.v_oq_r;
                

                if Control_type ==1
                    % AC voltage control
                    error_vd = v_dr - v_d- (i_gd*Rov-i_gq*Xov) * (-1);
              	    error_vq = v_qr - v_q- (i_gq*Rov+i_gd*Xov) * (-1);
                    i_d_r = -(error_vd*kpv1 + v_d_i);
                    i_q_r = -(error_vq*kpv1 + v_q_i);
                    dv_d_i = error_vd*kiv1;
                    dv_q_i = error_vq*kiv1;
    
                    % AC current control
                    error_id = i_d_r-i_d;
                    error_iq = i_q_r-i_q;
                    e_d = -(error_id*kpi1 + i_d_i);
                    e_q = -(error_iq*kpi1 + i_q_i);
                    di_d_i = error_id*kii1;            
                    di_q_i = error_iq*kii1;
                else
                    dv_d_i = 0;
                    dv_q_i = 0;
                    di_d_i = 0;
                    di_q_i = 0;
                    e_d   =v_dr;
                    e_q   =v_qr;
                end

                Pr = i_l * ed_dc;
                p  = (e_d*i_d + e_q*i_q)*(-1);
                w = ((Pr - p)/v_o*R1 + (v_o-V_dc))*N1*W0+w0; 

                % dtheta = w-W0;
                dtheta = w;
                dv_o = (Pr-p)/v_o/Co; 

                % Lf equation
                % e_d - v_od = -(di_ld/dt*Lf + Rf*i_ld - w*Lf*i_lq)
                % e_q - v_oq = -(di_lq/dt*Lf + Rf*i_lq + w*Lf*i_ld)
                di_d = (v_d - e_d - Rf*i_d + w*Lf*i_q)/Lf;
                di_q = (v_q - e_q - Rf*i_q - w*Lf*i_d)/Lf;

                % Cf equation
                % -(i_ld - i_od) = Cf*dv_cd/dt - w*Cf*v_cq
                % -(i_lq - i_oq) = Cf*dv_cq/dt + w*Cf*v_cd
                dv_d = (-(i_d - i_gd) + w*Cf*v_q)/Cf;
                dv_q = (-(i_q - i_gq) - w*Cf*v_d)/Cf;

                % Lc equation
                % v_od - v_d = -(Lc*di_od/dt + Rc*i_od - w*Lc*i_oq)
                % v_oq - v_q = -(Lc*di_oq/dt + Rc*i_oq + w*Lc*i_od)
                di_gd = (vgd - v_d - Rc*i_gd + w*Lc*i_gq)/Lc;
                di_gq = (vgq - v_q - Rc*i_gq - w*Lc*i_gd)/Lc;
                
                % dx
                f_xu = [dv_i; di_i; di_l; dv_pv; dv_d_i; dv_q_i; di_d_i; di_q_i; di_d; di_q; dv_d; dv_q; di_gd; di_gq; dv_o; dtheta];
                Output = f_xu;
                
            elseif CallFlag == 2    
                error_v_pv = vr-v_pv;
                i_r        = kpv_dc*error_v_pv + v_i;

                error_i_l  = i_r - i_l;
                ed_dc        = kpi_dc*error_i_l + i_i;
                
                v_dr = obj.v_od_r;
                v_qr = obj.v_oq_r;
                
                if Control_type ==1
                    % AC voltage control
                    error_vd = v_dr - v_d- (i_gd*Rov-i_gq*Xov) * (-1);
              	    error_vq = v_qr - v_q- (i_gq*Rov+i_gd*Xov) * (-1);
                    i_d_r = -(error_vd*kpv1 + v_d_i);
                    i_q_r = -(error_vq*kpv1 + v_q_i);
                    
                    % AC current control
                    error_id = i_d_r-i_d;
                    error_iq = i_q_r-i_q;
                    e_d = -(error_id*kpi1 + i_d_i);
                    e_q = -(error_iq*kpi1 + i_q_i);
                else
                    e_d   = v_dr;
                    e_q   = v_qr;
                end

                % end
                Pr = i_l * ed_dc;
                p  = (e_d*i_d + e_q*i_q)*(-1);
                w = (v_o+ (Pr - p)/v_o*R1 - 1)*N1*W0+w0; 

                % dx
                g_xu = [i_gd; i_gq; w; theta];
                Output = g_xu;
            end
            
        end
        
    end
end




