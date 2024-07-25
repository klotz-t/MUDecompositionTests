classdef MotorNeuronPool_CisiKohn_2008 
    % This is an implementation of the Motor Neuron Model of Cisi & Kohn
    % (2008). Thereby, the synaptic conductances are not considered in this
    % implementation. Further, the pulse-based approximations for the ion
    % channels have been removed and are rather represented by the full
    % governing differential equations. The gating dynamics of the ion
    % channels are modelled as in Negro & Farina (2011), based on Traub et
    % al (1991).
    % 
    % Literature:
    % - Cisi & Kohn (2008): Simulation system of spinal cord motor nuclei and
    %       associated nerves and muscles, in a Web-based architecture.
    % - Negro & Farina (2011): Decorrelation of cortical inputs and motoneuron output.
    % - Traub, Wong, Miles & Michelson (1991): A model of a CA3 hippocampal
    %       pyramidal neuron incorporating voltage-clamp data on intrinsic conductances.
    % - Heidlauf (2015):  Chemo-electro-mechanical modelling of the neuromuscular system.
    
    properties
        NumberOfMUs 
        Parameters
        InitStates
        Input
    end
    
    methods
        % -------------------------------------------------------------- %
        % -------- Evaluate a single motorneuron from the pool --------- %
        function dy = evaluate_rhs_single_mn(obj,time,y,mu_idx,drive)
            % Evaluates the right hand side of the motor neuron model of
            % Cisi and Kohn 2008. Therein y is the state vector, mu_idx is
            % the Number of the Motor Neuron and drive represents all
            % Inputs of the model.
            % y(2) is membrane potential in mV and t time in ms
            dy(1) = (-obj.Parameters(1,mu_idx)*(y(1)-obj.Parameters(11,mu_idx))- ...
                obj.Parameters(5,mu_idx)*(y(1)-y(2)))/obj.Parameters(7,mu_idx);
            dy(2) = (drive-obj.Parameters(6,mu_idx)*(y(2)-obj.Parameters(11,mu_idx))-...
                obj.Parameters(5,mu_idx)*(y(2)-y(1))-obj.Parameters(4,mu_idx)*y(3).^3*y(4).*(y(2)-...
                obj.Parameters(9,mu_idx))-obj.Parameters(2,mu_idx)*y(5).^4*(y(2)-obj.Parameters(10,mu_idx))-...
                obj.Parameters(3,mu_idx)*y(6).^2*(y(2)-obj.Parameters(10,mu_idx)))/obj.Parameters(8,mu_idx);
            % Gating Variables
            dy(3) = 0.32*(13-y(2))/(exp((13-y(2))/5)-1)*(1-y(3))-0.28*(y(2)-40)/(exp((y(2)-40)/5)-1)*(y(3)); % % m-gate (Na+ channel) -- equation from Negro & Farina (2011) / Traub et al (1991)
            dy(4) = 0.128*(exp((17-y(2))/18))*(1-y(4))-4/(exp((40-y(2))/5)+1)*(y(4)); % h-gate (Na+ channel) -- equation from Negro & Farina (2011) 
            dy(5) = 0.032*(15-y(2))/(exp((15-y(2))/5)-1)*(1-y(5))-0.5*(exp((10-y(2))/40))*(y(5)); % n-gate (fast K+ channel) -- equation from Negro & Farina (2011) / Traub et al (1991)
            dy(6) = 3.5/(exp((55-y(2))/4)+1)*(1-y(6))-0.025*(y(6)); % q-gate (slow K+ channel) -- equation from Negro & Farina (2011)
            dy = dy';            
        end
        % -------------------------------------------------------------- %
        % --------- Evaluate a all motorneurons from the pool using a for loop ---------- %
%         function dy = evaluate_rhs_cisi_kohn_2008_all_MNs(obj,time,y,drive)
%             % Evaluates the right hand side of the motor neuron model of
%             % Cisi and Kohn 2008. Therein y is the state vector and drive
%             represents all % Inputs of the model. % y(2) is membrane
%             potential in mV and t time in ms
%             for mu_idx = 1:obj.NumberOfMUs
%                 start_int = (mu_idx-1)*6;
%                 dy(start_int+1) = (-obj.Parameters(1,mu_idx)*(y(start_int+1)-obj.Parameters(11,mu_idx))- ...
%                     obj.Parameters(5,mu_idx)*(y(start_int+1)-y(start_int+2)))/obj.Parameters(7,mu_idx);
%                 dy(start_int+2) = (drive-obj.Parameters(6,mu_idx)*(y(start_int+2)-obj.Parameters(11,mu_idx))-...
%                     obj.Parameters(5,mu_idx)*(y(start_int+2)-y(start_int+1))-obj.Parameters(4,mu_idx)*y(start_int+3).^3*y(start_int+4)*(y(start_int+2)-...
%                     obj.Parameters(9,mu_idx))-obj.Parameters(2,mu_idx)*y(start_int+5).^4*(y(start_int+2)-obj.Parameters(10,mu_idx))-...
%                     obj.Parameters(3,mu_idx)*y(start_int+6).^2*(y(start_int+2)-obj.Parameters(10,mu_idx)))/obj.Parameters(8,mu_idx);
%                 % Gating Variables
%                 dy(start_int+3) = 0.32*(13-y(start_int+2))/(exp((13-y(start_int+2))/5)-1)*(1-y(start_int+3))-...
%                     0.28*(y(start_int+2)-40)/(exp((y(start_int+2)-40)/5)-1)*(y(start_int+3)); % m
%                 dy(start_int+4) = 0.128*(exp((17-y(start_int+2))/18))*(1-y(start_int+4))-4/(exp((40-y(start_int+2))/5)+1)*(y(start_int+4)); % h
%                 dy(start_int+5) = 0.032*(15-y(start_int+2))/(exp((15-y(start_int+2))/5)-1)*(1-y(start_int+5))-...
%                     0.5*(exp((10-y(start_int+2))/40))*(y(start_int+5)); % n
%                 dy(start_int+6) = 3.5/(exp((55-y(start_int+2))/4)+1)*(1-y(start_int+6))-0.025*(y(start_int+6)); %q
%             end
%             dy = dy';            
%         end
        % -------------------------------------------------------------- %
        % --------- Evaluate a all motorneurons from the pool without for loop ---------- %
        function dy = evaluate_rhs_mn_pool(obj,time,y,drive)
            % Evaluates the right hand side of the motor neuron model of
            % Cisi and Kohn 2008. Therein y is the state vector and drive
            % represents all Inputs of the model.
            % y(2) is membrane potential in mV and t time in ms;
            y=reshape(y,6,[]);
                dy(1,:) = (-obj.Parameters(1,:).*(y(1,:)-obj.Parameters(11,:))- ...
                    obj.Parameters(5,:).*(y(1,:)-y(2,:)))./obj.Parameters(7,:);
                dy(2,:) = (drive-obj.Parameters(6,:).*(y(2,:)-obj.Parameters(11,:))-...
                    obj.Parameters(5,:).*(y(2,:)-y(1,:))-obj.Parameters(4,:).*y(3,:).^3.*y(4,:).*(y(2,:)-...
                    obj.Parameters(9,:))-obj.Parameters(2,:).*y(5,:).^4.*(y(2,:)-obj.Parameters(10,:))-...
                    obj.Parameters(3,:).*y(6,:).^2.*(y(2,:)-obj.Parameters(10,:)))./(obj.Parameters(8,:));
                % Gating Variables
                dy(3,:) = 0.32.*(13-y(2,:))./(exp((13-y(2,:))./5)-1).*(1-y(3,:))-...
                    0.28.*(y(2,:)-40)./(exp((y(2,:)-40)./5)-1).*(y(3,:)); % m-gate (Na+ channel) -- equation from Negro & Farina (2011) / Traub et al (1991) 
                dy(4,:) = 0.128.*(exp((17-y(2,:))./18)).*(1-y(4,:))-4./(exp((40-y(2,:))./5)+1).*(y(4,:)); % h-gate (Na+ channel) -- equation from Negro & Farina (2011) 
                %dy(4,:) = 0.06.*exp((-55-y(2,:)+75)/15).*(1-y(4,:))-6.01.*(exp((17-y(2,:)+75)/21)+1).*y(4,:);
                dy(5,:) = 0.032.*(15-y(2,:))./(exp((15-y(2,:))./5)-1).*(1-y(5,:))-...
                    0.5.*(exp((10-y(2,:))./40)).*(y(5,:)); % n-gate (fast K+ channel) -- equation from Negro & Farina (2011) / Traub et al (1991)
             dy(6,:) = 3.5./(exp((55-y(2,:))./4)+1).*(1-y(6,:))-0.025.*(y(6,:)); % q-gate (slow K+ channel) -- equation from Negro & Farina (2011) 
            %     dy(6,:) = 0.6.*3.5./(exp((47-y(2,:))./4)+1).*(1-y(6,:))-0.6.*0.025.*(y(6,:)); % g-gate tests for alternatives   
             % dy(6,:) = (1./(exp((y(2,:)+31-75)./(-15))+1)-y(6,:))./(5./(exp((y(2,:)+50-75)./40)+exp(-(y(2,:)+50-75)./50))); % ElBasiouny 2005
                % dy(6,:) = (1./(1+exp(-(y(2,:)+25-75)./20))-y(6,:))./(0.8+20*exp((y(2,:)+39-75)/5.5)./(1+exp((y(2,:)+39-75)./5.5)).^2); % formulation by Powers et al. 2012
%   dy(7,:) = 0.001.*exp((-85-y(2,:)+75)/30).*(1-y(7,:))-0.0034.*(exp((-17-y(2,:)+75)/10)+1).*y(7,:);
                dy = reshape(dy,[],1); 
        end
        % -------------------------------------------------------------- %
        function obj = init_motor_neuron_pool(obj,N)
            % Set Number of Motor Neurons and initialise all model parameters
            % by default values.
            obj.NumberOfMUs = N;
            % Initialise Model Parameters and define Initial Conditions [Cisi and Kohn 2008]
            if(obj.NumberOfMUs > 0)
                % Rows specify a specific material parameter
                % Columns represent the properties of a specific motor
                % neuron.
                obj.Parameters = zeros(11,N);
                obj.Parameters(1,:) = 5.1194e-04;         % Gld (1/kOhm)
                obj.Parameters(2,:) = 7.6181e-04;         % Gkf (1/kOhm) 
                obj.Parameters(3,:) = 0.003;              % Gks (1/kOhm) 
                obj.Parameters(4,:) = 0.0057;             % Gna (1/kOhm)  
                obj.Parameters(5,:) = 7.1070e-04;         % Gc (1/kOhm)
                obj.Parameters(6,:) = 1.6634e-04;         % Gls (1/kOhm)
                obj.Parameters(7,:) = 0.0073;             % Cd (uF)
                obj.Parameters(8,:) = 1.9045e-04;         % Cs (uF) 
                obj.Parameters(9,:) = 120;                % V_Na (mV)
                obj.Parameters(10,:) = -10;               % E_K  (mV)
                obj.Parameters(11,:) = 0;                 % E_L  (mV)
                % Rows specify a specific state variable, Columns represent 
                % the properties of a specific motor neuron.
                obj.InitStates = zeros(6,N);
                obj.InitStates(1,:) = 0.0027;             % Vm_d
                obj.InitStates(2,:) = 0.0579;             % Vm_s  
                obj.InitStates(3,:) = 0.0293;             % m  
                obj.InitStates(4,:) = 0.9959;             % h  
                obj.InitStates(5,:) = 0.0381;             % n
                obj.InitStates(6,:) = 1.5155e-04;         % q
               % obj.InitStates(7,:) = 0.0034;
            else
               error('Number of MUs is not defined correctly') 
            end
        end
        % -------------------------------------------------------------- %
        function obj = set_model_parameter(obj,value,para_idx,mu_idx)
            % Set a specific model parameter, defined by para_idx and 
            % mu_idx to the specified value.
            obj.Parameters(para_idx,mu_idx) = value; 
        end
        % -------------------------------------------------------------- %
        function obj = set_init_state(obj,value,state_idx,mu_idx)
            % Set a specific inital state, defined by state_idx and 
            % mu_idx to the specified value.
            obj.Parameters(state_idx,mu_idx) = value; 
        end
        % -------------------------------------------------------------- %
        function obj = distribute_parameters(obj,Cm,Ri,Rmd,Rms,...
                                                        ld,ls,rd,rs,type,a)
            % Distribute the model parameters between two specified extreme
            % values. As default setting it is assumed that the 
            % Parameters Rmd,Rms,ld,ls,rd and rs are exponetially 
            % distributed (cf. Negro & Farina (2011). 
            % Further, Cm and Ri are constant for all MNs.
            N = obj.NumberOfMUs;
            
            if nargin == 9 
                type = 'exp'; % Default: exponential distribution
                spec = 100;   % Parameter specifying the shape of the exp distribution
            elseif nargin == 11
                spec = a;
            end
            
            if (strcmp(type,'exp'))
                Rmd_values = obj.exp_distribution(Rmd(1),Rmd(2),spec); % Rmd (kOhm cmÂ²)
                Rms_values = obj.exp_distribution(Rms(1),Rms(2),spec); % Rms (kOhm cmÂ²)
                ld_values  = obj.exp_distribution(ld(1),ld(2),spec); % ld (cm)
                ls_values  = obj.exp_distribution(ls(1),ls(2),spec); % ls (cm)
                rd_values  = obj.exp_distribution(rd(1),rd(2),spec); % rd (cm)
                rs_values  = obj.exp_distribution(rs(1),rs(2),spec); % rs (cm)
            elseif (strcmp(type,'lin'))
                Rmd_values = linspace(Rmd(1),Rmd(2),N);   % Rmd (kOhm cmÂ²)
                Rms_values = linspace(Rms(1),Rms(2),N);   % Rms (kOhm cmÂ²)
                ld_values  = linspace(ld(1),ld(2),N);     % ld (cm)
                ls_values  = linspace(ls(1),ls(2),N);     % ls (cm)
                rd_values  = linspace(rd(1),rd(2),N);     % rd (cm)
                rs_values  = linspace(rs(1),rs(2),N); % rs (cm)
            elseif (strcmp(type,'uniform'))
                Rmd_values = Rmd.*ones(1,N);   % Rmd (kOhm cmÂ²)
                Rms_values = Rms.*ones(1,N);   % Rms (kOhm cmÂ²)
                ld_values  = ld.*ones(1,N);     % ld (cm)
                ls_values  = ls.*ones(1,N);     % ls (cm)
                rd_values  = rd.*ones(1,N);     % rd (cm)
                rs_values  = rs.*ones(1,N); % rs (cm)
            else
                fprintf('Error: Specified distribution type is not implemented')
            end
            
            % Store geometric parameters in object
            obj.Input(1,:) = Cm.*ones(1,N);
            obj.Input(2,:) = Ri.*ones(1,N);
            obj.Input(3,:) = Rmd_values;
            obj.Input(4,:) = Rms_values;
            obj.Input(5,:) = ld_values;
            obj.Input(6,:) = ls_values;
            obj.Input(7,:) = rd_values;
            obj.Input(8,:) = rs_values;
            
            obj.Parameters(1,:) = 2*pi.*rd_values.*ld_values./Rmd_values;   % Gld (1/kOhm)
            obj.Parameters(2,:) = 4*2*pi.*rs_values.*ls_values;             % Gkf (1/kOhm) (max. conductance = 4 mS/cmÂ²)
            obj.Parameters(3,:) = 16*2*pi.*rs_values.*ls_values;            % Gks (1/kOhm) (max. conductance = 16 mS/cmÂ²)
            obj.Parameters(4,:) = 30*2*pi.*rs_values.*ls_values;            % Gna (1/kOhm) (max. conductance = 30 mS/cmÂ²)
            obj.Parameters(5,:) = 2./(Ri.*ld_values./(pi.*rd_values.^2)+ ...
                Ri.*ls_values./(pi.*rs_values.^2));                         % Gc (1/kOhm)
            obj.Parameters(6,:) = 2*pi.*rs_values.*ls_values./Rms_values;   % Gls (1/kOhm)
            obj.Parameters(7,:) = 2*pi.*rd_values.*ld_values.*Cm;           % Cd (uF)
            obj.Parameters(8,:) = 2*pi.*rs_values.*ls_values.*Cm;           % Cs (uF)
            obj.Parameters(9,:) = 120;                                      % V_Na (mV)
            obj.Parameters(10,:) = -10;                                     % E_K  (mV)
            obj.Parameters(11,:) = 0;                                       % E_L  (mV)            
        end
        % -------------------------------------------------------------- %
        function vector = exp_distribution(obj,value1,value2,p1)
            % Generates an exponentially distributed vector of dimension N 
            % in the intervall specified by 'value1' and 'value2'.
            N = obj.NumberOfMUs;
            intervall = value2 -value1;
            gg = linspace(0,1,N);   % Distribute MN indexes between 0 and 1
            vector = exp(log(p1).*gg)./p1.*intervall+value1;
        end
          

    end
    
end



