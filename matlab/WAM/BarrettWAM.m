% Anthropomorphic arm (3shoulder+1elbow+1wrist)

classdef BarrettWAM < handle

    properties   
        % disturbance
        dist
        % use cpp for dynamics models
        cpp
        % actual parameter structure
        PAR_ACT
        % nominal parameters structure
        PAR
        % constraints structure
        CON
        % fields necessary for simulation and plotting, noise etc.
        SIM       
        % forward nominal and actual dynamics handles
        % used if cpp is true
        f_forward_nom
        f_forward_act
        % inverse dynamics handle used if cpp is true
        f_inverse
    end
    
    methods
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            
            % check that the input has all the fields
            %assert(all(strcmp(fieldnames(obj.CON), fieldnames(STR))));
            obj.CON = STR;
            
        end 
        
        % set the simulation parameters
        function set.SIM(obj, sim)
            obj.SIM.discrete = sim.discrete;
            obj.SIM.dimx = sim.dimx;
            obj.SIM.dimy = sim.dimy;
            obj.SIM.dimu = sim.dimu;
            obj.SIM.h = sim.h;
            obj.SIM.eps_m = sim.eps_m;
            %assert(strcmpi(sim.int,'Euler') || strcmpi(sim.int,'RK4'),...
            %       'Please input Euler or RK4 as integration method');
            obj.SIM.int = sim.int;
            obj.SIM.C = sim.C;
        end
        
        % disturbance function
        % uncomment lines in actual() to use different wam link parameters
        function set_dist(obj)
           
            Nmax = 1000;
            dt = obj.SIM.h;
            t = dt * (1:Nmax);
            obj.dist =  [zeros(7,Nmax); 5.0 * sample_traj(7,t,0.1,zeros(7,1))']; 
        end
        
    end
    
    methods
        
        % constructor for convenience
        function obj = BarrettWAM(con,sim,dyn)
            
            
            obj.cpp = dyn.use_cpp;
            obj.SIM = sim;
            % set object parameter
            obj.PAR =  load_wam_links(dyn.nom);
            obj.PAR_ACT = load_wam_links(dyn.act);
            % set object constraints
            obj.CON = con;
            % repeating/not repeating disturbances
            %obj.set_dist();
            
            link_nom = dyn.nom;
            link_act = dyn.act;
            obj.f_forward_nom = @(q,qd,u) barrett_wam_dyn_art(link_nom,q,qd,u);
            obj.f_inverse = @(q,qd,qdd) barrett_wam_invdyn_ne(link_nom,q,qd,qdd);
            obj.f_forward_act = @(q,qd,u) barrett_wam_dyn_art(link_act,q,qd,u);
        end
        
        %% Dynamics functions here
        % provides nominal model
        function x_dot = nominal(obj,~,x,u)
            % differential equation of the forward dynamics
            % x_dot = A(x)x + B(x)u + C(x)
            if obj.cpp
                x_dot = [x(8:end);obj.f_forward_nom(x(1:7),x(8:end),u)];
            else
                x_dot = barrettWamDynamicsArt(x,u,obj.PAR);
            end
        end
        
        % provides actual model
        function x_dot = actual(obj,t,x,u)
            
            % differential equation of the inverse dynamics
            % x_dot = A(x) + B(x)u
            dt = obj.SIM.h;
            if obj.cpp                
                %x_dot = [x(8:end);barrett_wam_dyn_art_nom(x(1:7),x(8:end),u)];
                %x_dot = x_dot + obj.dist(:,ceil(t/dt));
                x_dot = [x(8:end);obj.f_forward_act(x(1:7),x(8:end),u)];
            else
                x_dot = barrettWamDynamicsArt(x,u,obj.PAR_ACT);            
                %x_dot = barrettWamDynamicsArt(x,u,obj.PAR) + obj.dist(:,ceil(t/dt));            
            end
            
        end
        
        % dynamics to get u
        % TODO: can we input/output matrices?
        function u = invDynamics(obj,q,qd,qdd)
            % inverse dynamics model taken from SL
            if obj.cpp
                u = obj.f_inverse(q,qd,qdd);
            else
                u = barrettWamInvDynamicsNE(q,qd,qdd,obj.PAR);
                %u = barrettWamInvDynamicsArt(q,qd,qdd,obj.PAR);
            end
        end
        
        % linearizes the nominal dynamics around the trajectory
        function [A,B,Ac,Bc] = linearize(obj,s,uff)
            
            N = size(uff,2);
            h = obj.SIM.h;
            %assert(size(uff,2) == N, 'controls and traj should have same length!');
            dimx = obj.SIM.dimx;
            dimu = obj.SIM.dimu;
            %h = obj.SIM.h;
            A = zeros(dimx,dimx,N);
            B = zeros(dimx,dimu,N);
            Ac = zeros(dimx/2,dimx,N); % qdd = A[q;qd] + Bu;
            Bc = zeros(dimx/2,dimu,N);
            for i = 1:N
                if obj.cpp
                    [Amat,Bmat] = num_diff(obj.f_forward_nom,s(:,i),uff(:,i),obj.PAR,true);
                else
                    [Amat,Bmat] = num_diff(@barrettWamDynamicsArt,s(:,i),uff(:,i),obj.PAR,false);
                end
                % get discrete approximation from jacobian
                % crude approximation
                %A(:,:,i) = eye(dimx,dimx) + obj.SIM.h * A(:,:,i);
                %B(:,:,i) = obj.SIM.h * B(:,:,i);
                % exact matrix calculation 
                %%{
                Mat = [Amat, Bmat; zeros(dimu, dimx + dimu)];
                %Mat = [A(:,:,i), B(:,:,i); zeros(dimu, dimx + dimu)];
                MD = expm(h * Mat);
                A(:,:,i) = MD(1:dimx,1:dimx);
                B(:,:,i) = MD(1:dimx,dimx+1:end);
                
                Ac(:,:,i) = Amat(dimx/2+1:end,:);
                Bc(:,:,i) = Bmat(dimx/2+1:end,:);
                %}
            end
        end
        
        %% Kinematics related functions here
        % Calculate the racket orientation based on quaternion
        function racketOrient = calcRacketOrientation(obj,cartOrient)
            
            % quaternion transformation of -pi/2 from endeff to racket
            rot = [cos(pi/4); -sin(pi/4); 0; 0];
            
            racketOrient = mult2Quat(cartOrient,rot);   
        end
        
        % run kinematics using an external function
        % return endeffector coordinates, vel and orientation
        % TODO: should return what barrettWAMKinematics returns
        % make sure endeffector[Z] = 0
        function [x,xd,o] = kinematics(obj,Q)
              
            assert(size(Q,1) == 14, 'velocities not fed in!');
            dim = size(Q,1)/2;
            lenq = size(Q,2);
            q = Q(1:dim,:);
            qd = Q(dim+1:end,:);
            x = zeros(3,lenq);
            xd = zeros(3,lenq);
            o = zeros(4,lenq);
            for i = 1:lenq
                [xLink,xOrigin,xAxis,Amats] = barrettWamKinematics(q(:,i),obj.PAR);
                o(:,i) = rot2Quat(Amats(6,1:3,1:3));
                jac = jacobian(xLink,xOrigin,xAxis);
                x(:,i) = xLink(6,1:3)';
                xd(:,i) = jac(1:3,:) * qd(:,i);
            end
        end
        
        % racket configurations and velocities are returned
        % racket orientations are also returned as a quaternions        
        function [x,xd,o] = calcRacketState(obj,q,qd)
            
            lenq = size(q,2);
            x = zeros(3,lenq);
            xd = zeros(3,lenq);
            o = zeros(4,lenq);
            for i = 1:lenq
                [xLink,xOrigin,xAxis,Amats] = barrettWamKinematics(q(:,i),obj.PAR);
                quat = rot2Quat(Amats(6,1:3,1:3));
                o(:,i) = obj.calcRacketOrientation(quat);
                jac = jacobian(xLink,xOrigin,xAxis);
                x(:,i) = xLink(6,1:3)';
                xd(:,i) = jac(1:3,:) * qd(:,i);
            end
        end
        
        % calculate racket normal from racket orientations
        % NOTE:
        % what can be done directly is the following rotation
        % i.e. directly, without using quaternions
        % R1 = squeeze(Amats(6,1:3,1:3));
        % R2 = [1 0 0 ; 0 0 1; 0 -1 0]; 
        % %since theta = -pi/2, using
        % %angle/axis transformation [see Siciliano Robotics book pg.54]
        % R = R1 * R2;
        % normal = R(:,3);
        function normal = calcRacketNormal(obj,o)
            rot = quat2Rot(o);
            normal = rot(1:3,3);
        end
        
        % calculate the geometric jacobian using the common jacobian function
        function jac = calcJacobian(obj,q)
            
            [xLink,xOrigin,xAxis,~] = barrettWamKinematics(q,obj.PAR);
            jac = jacobian(xLink,xOrigin,xAxis);
        end
        
        %% Inverse Kinematics functions
        
        % Inverse Kinematics for table tennis
        % Specifically, trying to keep initial slide of the racket
        function [qf,qfdot] = invKinTableTennis(obj,Q0,racket)
            
            racketPos = racket.pos;
            racketNormal = racket.normal;
            racketVel = racket.vel;
            racketAngularVel = racket.angvel;
            q0 = Q0(1:7);  
            q0d = Q0(8:end);
            
            try
                tic;
                % get the slide of the original racket orientation
                [~,~,o] = obj.calcRacketState(q0,q0d);
                rotMatrix0 = quat2Rot(o);
                slide0 = rotMatrix0(1:3,2);    
                slide0 = slide0./norm(slide0,2);
                % add a comfortable slide close to original slide
                a = racketNormal;
                % project slide0 to the racket plane
                %projNormal = racketNormal*racketNormal';
                %projRacket = eye(3) - projNormal;
                %s = projRacket * slide0;    
                s = slide0 - (racketNormal'*slide0).*racketNormal;    
                %assert(abs(norm(s,2)-1) < 1e-2,'slide not normalized');
                s = s./norm(s,2);
                assert(s'*a < 1e-3,'slide calculation not working!');
                n = crossProd(s,a);                
                rotMatrix = [n,s,a];
                quatRacket = rot2Quat(rotMatrix);
                % get endeffector position and quaternion
                ePos = racketPos;
                rotBack = [cos(-pi/4); -sin(-pi/4); 0; 0];
                eQuat = mult2Quat(quatRacket,rotBack);
                % estimate qf from demonstrations
                %qest = [ePos(:);eQuat(:)]' * obj.Bdemo;
                %qest = obj.clampJointLimits(qest(:));
                qest = rand(7,1);
                qf = obj.invKinematics(ePos(:),eQuat(:),qest(:));
                obj.checkJointLimits(qf,zeros(7,1),zeros(7,1));
                timeInvKin = toc;
                fprintf('InvKin took %f sec.\n',timeInvKin);            
                qfdot = obj.calcJacobian(qf) \ [racketVel;racketAngularVel];
            catch ME
                disp(ME.message);
                disp('InvKin problem. Not moving the robot...');                
                %disp('Virtual Hitting Point outside of workspace');
                qf = q0;
                qfdot = zeros(7,1);
            end
        end
        
        % call inverse kinematics from outside
        % o is a quaternion
        function q = invKinematics(obj,x,o,q0)
            
            % not masking any cartesian space variables
            %M = [1 1 1 1 1 1];
    
            fprintf('Computing IK for WAM...\n');
            if size(x,2) == 1
                % solve with fmincon
                q = invKinematicsOpt(obj,x,o,q0);
            else
                for k = 1:1:size(x,2) 

                    T = [quat2Rot(o(:,k)),x(:,k);zeros(1,3),1];
                    q(:,k) = invKinematicsGen(T, q0, obj.PAR);
                    q0 = q(:,k); % use previous joint values as initial guess for ikine
                end
            end
            %}
           
            
            %q = unwrap(q);
            
        end        

        %% Check safety here
        
        % TODO: check torques in a more clever way!
        function [q,qd,qdd] = checkJointLimits(obj,q,qd,qdd)
            
            len = size(q,2);
            con = obj.CON;
            umax = repmat(con.u.max,1,len);
            umin = repmat(con.u.min,1,len);
            qmax = repmat(con.q.max,1,len);            
            qmin = repmat(con.q.min,1,len);
            qdmax = repmat(con.qd.max,1,len);
            qdmin = repmat(con.qd.min,1,len);
            qddmax = repmat(con.qdd.max,1,len);
            qddmin = repmat(con.qdd.min,1,len); 
            %u = zeros(7,len);
            %for i = 1:len
            %    u(:,i) = obj.invDynamics(q(:,i),qd(:,i),qdd(:,i));
            %end
            
            %obj.displayMaxInfo(q,qd,qdd,u);
            
            try              
                assert(~any(any((q > qmax | q < qmin))),'Joint limits violated!');
                assert(~any(any((qd > qdmax | qd < qdmin))), 'Vel limits violated!');
                assert(~any(any((qdd > qddmax | qdd < qddmin))), 'Acc limits violated!');
                %assert(~any(any((u > umax | u < umin))), 'Torque limits violated!');
            catch ME
                disp(ME.message);
                disp('Not moving the robot!');
                dof = length(con.q.max);
                q0 = q(:,end);
                q = q0 * ones(1,len);
                qd = zeros(dof,len);
                qdd = zeros(dof,len);
            end
        end
        
        % Display max and min info about trajectory and control inputs
        function displayMaxInfo(obj,q,qd,qdd,u)
            
            qmax = max(q,[],2);
            qmin = min(q,[],2);
            qdmax = max(qd,[],2);
            qdmin = min(qd,[],2);
            qddmax = max(qdd,[],2);
            qddmin = min(qdd,[],2);
            umax = max(u,[],2);
            umin = min(u,[],2);
            
            fprintf('qmax 	 qdmax 	 qddmax    umax\n');
            for i = 1:7
                fprintf('%-8.2f %-8.2f %-8.2f %-8.2f\n',...
                    qmax(i),qdmax(i),qddmax(i),umax(i));
            end
            fprintf('qmin 	 qdmin 	 qddmin    umin\n');
            for i = 1:7
                fprintf('%-8.2f %-8.2f %-8.2f %-8.2f\n',...
                    qmin(i),qdmin(i),qddmin(i),umin(i));
            end
                
        end
        
        function q = clampJointLimits(obj,q)
            
            con = obj.CON;
            q = min(con.q.max,q);
            q = max(con.q.min,q);
        end
        
        function q = checkContactWithTable(obj,q)
            % TODO:
        end
        
        %% Drawing functions here
        % to draw the robots joints and endeffector 
        % for one posture only
        function [joints,endeff,racket] = drawPosture(obj,q)
            
            [xLink,xOrigin,xAxis,Amats] = barrettWamKinematics(q,obj.PAR);
            quat = rot2Quat(Amats(6,1:3,1:3));
            orient = obj.calcRacketOrientation(quat);
            R = quat2Rot(orient(:));
            joints = xOrigin;
            endeff = xLink(6,1:3);
            
            %x and y are the coordinates of the center of the circle
            %r is the radius of the circle
            %0.01 is the angle step, bigger values will draw the circle faster but
            %you might notice imperfections (not very smooth)
            racket_radius = 0.076;
            ang = 0:0.01:2*pi; 
            racket_x = racket_radius * cos(ang);
            racket_y = racket_radius * sin(ang);
            % transform into base coord.
            racket = repmat(endeff(:),1,length(ang)) + ...
                R * [racket_x; racket_y; zeros(1,length(ang))];
            
        end
        
        % make an animation of the robot manipulator
        function animateArm(obj,qs)
            
            %TODO:
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift_constraints(obj,trj,ilc)
            
            %TODO:
            
        end
        
        %% PLOTTING FUNCTIONS
        % plot the control inputs
        function plot_inputs(obj,U)
            
            N = size(U,2);
            dt = obj.SIM.h;
            t = dt * (1:N);
            num_iter = size(U,3);
            num_inp = size(U,1);
            figure;
            for k = 1:num_iter
                for i = 1:num_inp
                    subplot(num_inp,1,i);
                    plot(t,U(i,:,k));
                    %title(strcat(num2str(i),'. control input'));
                    ylabel(['$u_',num2str(i),'$'],'Interpreter','latex');
                end
                xlabel('Time (s)');
            end
            
        end
        
        % show experiment outcomes
        % experiment means that we average ILC learning stats
        % over different trajectories
        function show_experiment(obj,exp)
            
            % experiments should all have same num of iterations
            num_exp = length(exp);
            num_iter = size(exp{1}.x,3);
            err_norm = zeros(1,num_iter);
            final_cost = zeros(num_iter,num_exp);            
            labels_x = cell(1,4);
            labels_u = cell(1,3);

            for i = 1:num_iter
                for j = 1:num_exp
                    e = exp{j}.x(:,:,i) - exp{j}.r(:,:);
                    err_norm(i,j) = norm(e,'fro');
                    final_cost(i,j) = norm(e(:,end));
                end
                labels_u{i} = ['trial ', int2str(i)];
                labels_x{i+1} = ['trial ', int2str(i)];
            end
            labels_x{1} = 'desired traj';
            
            figure;
            subplot(2,1,1);
            if num_exp == 1                
                plot(0:1:num_iter-1,err_norm(:,1));
            else
                sigmas = std(err_norm,1,2);
                errorbar(mean(err_norm,2),sigmas);
            end
            ylabel('Error norm');
            xlabel('Iterations');
            title('Error norm vs trials');
            
            subplot(2,1,2);
            if num_exp == 1
                plot(1:num_iter,final_cost(:,1));
            else
                sigmas = std(final_cost,1,2);
                errorbar(mean(final_cost,2),sigmas);
            end
            ylabel('Final cost norm');
            xlabel('Iterations');
            title('Final cost norm vs trials');            
        end
        
        % plot the desired system states and outputs
        % if plot_all_cart is turned on all iterations will be
        % included in the cartesian plot
        function plot_outputs(obj,X,traj,plot_all_cart)
 
            N = size(X,2);
            dt = obj.SIM.h;
            t = dt * (1:N);
            num_iter = size(X,3);
            num_out = size(X,1);
            figure; % draw joint outputs
            
            for i = 1:num_out/2
                subplot(num_out/2,2,2*i-1);
                plot(t,X(i,:,end),'.-',t,traj(i,:),'--');
                if i == num_out/2
                    xlabel('Time (s)');
                end
                ylabel(['$q_', num2str(i), '$'], 'Interpreter', 'latex');
                subplot(num_out/2,2,2*i);
                plot(t,X(i+num_out/2,:,end),'.-',t,traj(i+num_out/2,:),'--');
                %legend(name,'Reference');
                %title(strcat(num2str(i),'. state'));
                %xlabel('Time (s)');
                ylabel(['$\dot{q}_', num2str(i), '$'], 'Interpreter', 'latex');
            end
            xlabel('Time (s)');
            hold on;            
            
            for k = 1:num_iter
                % display SSE error
                sse = norm(X(:,:,k) - traj, 'fro')^2;
                fincost = norm(X(:,end,k) - traj(:,end));
                rms = sqrt(sse/N);
                fprintf('RMS for ILC is %f \n', rms);
                fprintf('Final cost is %f \n', fincost);
            end

            figure; % draw Cartesian output
            yDesCart = obj.kinematics(traj);
            plot3(yDesCart(1,:),yDesCart(2,:),yDesCart(3,:),'r--');
            hold on;
            grid on;
            
            if plot_all_cart
                idx_draw = 1:num_iter;
            else
                idx_draw = num_iter;
            end
            
            for i = idx_draw
                yCart = obj.kinematics(X(:,:,i));
                plot3(yCart(1,:),yCart(2,:),yCart(3,:),'-');
            end
            xlabel('x');
            ylabel('y');
            zlabel('z');
            for i = 1:num_iter
                labels_x{i+1} = ['trial ', int2str(i)];
            end
            labels_x{1} = 'desired';
            legend(labels_x{1:4},'Location','northwest');
            %legend('desired','actual');

            % include time labels on Cart. traj.
            %%{
            drawTimeIter = 40;
            tLabel = t(1:drawTimeIter:end);
            precision = 4;
            tLabelCell = num2cell(tLabel,precision);
            for i = 1:length(tLabelCell)
                tLabelCell{i} = num2str(tLabelCell{i});
            end
            % annotate some of the ball positions
            xDraw = yDesCart(1,1:drawTimeIter:end);
            yDraw = yDesCart(2,1:drawTimeIter:end);
            zDraw = yDesCart(3,1:drawTimeIter:end);
            %text(xDraw,yDraw,zDraw,tLabelCell);
            scatter3(xDraw,yDraw,zDraw,20,'b','*');
            %}
            hold off;            
            
        end
        
    end
end