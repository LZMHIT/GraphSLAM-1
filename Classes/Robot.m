classdef Robot < handle & matlab.mixin.Copyable
    
    properties %(Access = private)
        iteration = 0;
        number = 1;
        numRobots = 1;
        states = [0;0;0];
        map = [];
        mapPointIdx = [];
        measurements = [];
        measCorrespond = [];
        lastTransmitIter = 0;
        anchor = [0;0;0];
        encounters = [];        
        control = [0;0];
        reorderRate = 100;
        dt = 0.001;
        Q = 0;
        R = 0;
        sqrtInfoR;
        rhsD;
        
        anchorIdx = [1:3]';
        stateIdx = [4:6]';
        mapIdx = [];
        
        
        hasPrior = false;
        priorMean = [];
        priorCov  = [];
        
        % Scaling factors for sqrt info matrix
        invSqrtQ;
        invSqrtR;
        invSqrtPriorCov;
        
        % Control traj
        controlTraj;
    end
    
    methods
        % Constructor
        function obj = Robot(numRobots,anchors,dt,Q,R,reorderRate,r1PriorMean,r1PriorCov,controlTraj,measurements,measCorrespond)          
            if nargin > 0
                obj = copy(repelem(obj,1,numRobots));
                
                if size(anchors,2) == 1
                    anchors = anchors.*ones(3,numRobots);
                end
                
                % Set prior for first robot
                obj(1,1).hasPrior = true;
                obj(1,1).priorMean = r1PriorMean;
                obj(1,1).priorCov = r1PriorCov;
                obj(1,1).invSqrtPriorCov = sqrtm(inv(r1PriorCov))';
                
                for i = 1:numRobots
                    obj(1,i).number = i;
                    obj(1,i).dt = dt;
                    obj(1,i).Q = Q;
                    obj(1,i).R = R;
                    obj(1,i).reorderRate = reorderRate;
                    obj(1,i).anchor = anchors(:,1);
                    
                    obj(1,i).invSqrtQ = sqrtm(inv(Q))';
                    obj(1,i).invSqrtR = sqrtm(inv(R))'; 
                    obj(1,i).controlTraj = controlTraj(:,:,i);
                    obj(1,i).measurements = measurements(:,:,i);
                    obj(1,i).measCorrespond = measCorrespond(:,:,i);
                end
            end            
        end
        
        % getters
        function c = getDeepCopy(obj)
            c = copy(obj);
        end
        
        % Motion model
        function xp = predictState(obj,x,u)
           xp = x + obj.dt*[u(1)*cos(x(3)); u(1)*sin(x(3)); u(2)]; 
        end
        % Measurement model
        function g = predictMeasurement(obj,x,m)
            g = [norm(m - x(1:2),2); % range
                 (atan2(m(2) - x(2),m(1) - x(1)) - x(3))]; % bearing 
        end
        
        % Get jacobians
        function F = getProcessJacobian(obj,x,u)
            F = [1  0  -u(1)*sin(x(3))*obj.dt;
                 0  1   u(1)*cos(x(3))*obj.dt;
                 0  0   1];
        end
        function [H,J] = getMeasurementJacobians(obj,x,m)
            % measurement type (1 = range, 2 = bearing)                      
            difference = m - x(1:2);
            rangeSq = sum(difference.^2);
            range = sqrt(rangeSq);
            
            % Jacobian wrt state
            H = [-difference(1)/range,  -difference(2)/range,    0;
                  difference(2)/rangeSq,-difference(1)/rangeSq, -1];
            
            % Jacobian wrt map            
            J = [ difference(1)/range,  difference(2)/range;
                 -difference(2)/rangeSq,difference(1)/rangeSq];
        end
        function [Cx1,Cx2,Cd1,Cd2] = getEncounterJacobians(obj,x1,x2,d1,d2)
            % Get pose sums
            poseSum1 = addPose(d1,x1);
            poseSum2 = addPose(d2,x2);
            
            % Get rotation matrices
            s1 = sin(d1(3));
            c1 = cos(d1(3));
            rot1 = [ c1 s1;
                    -s1 c1];            
            s2 = sin(d2(3));
            c2 = cos(d2(3));
            rot2 = [ c2 s2;
                    -s2 c2];
            
            % measurement type (1 = range, 2 = bearing)                      
            difference = poseSum2(1:2) - poseSum1(1:2);
            rangeSq = sum(difference.^2);
            range = sqrt(rangeSq);
            
            % Rotated diffs
            diffR1 = rot1*difference;
            diffR2 = rot2*difference;
            
            % Third column for Cd1/Cd1
            thirdColFactor1 = [x1(1)*s1 + x1(2)*c1; -x1(1)*c1 + x1(2)*s1];
            thirdColFactor2 = [-x2(1)*s2 - x2(2)*c2; x2(1)*c2 - x2(2)*s2];
            thirdCol1 = [(thirdColFactor1'*difference)./range;
                         (thirdColFactor1(2)*difference(1) - thirdColFactor1(1)*difference(2))/rangeSq - 1];
            thirdCol2 = [(thirdColFactor2'*difference)./range;
                         (thirdColFactor2(2)*difference(1) - thirdColFactor2(1)*difference(2))/rangeSq];
            
            % Jacobian wrt x1
            Cx1 = [-diffR1(1)/range,  -diffR1(2)/range,    0;
                    diffR1(2)/rangeSq,-diffR1(1)/rangeSq, -1];
                
            % Jacobian wrt x2
            Cx2 = [ diffR2(1)/range,   diffR2(2)/range,   0;
                   -diffR2(2)/rangeSq, diffR2(1)/rangeSq, 0];    
               
            % Jacobian wrt d1
            Cd1 = [-difference(1)/range,  -difference(2)/range,   thirdCol1(1);
                    difference(2)/rangeSq,-difference(1)/rangeSq, thirdCol1(2)];
                
            % Jacobian wrt d2
            Cd2 = [ difference(1)/range,   difference(2)/range,   thirdCol2(1);
                   -difference(2)/rangeSq, difference(1)/rangeSq, thirdCol2(2)];   
      
        end


        
        % Run one time step of functionality
        function i = runIteration(obj)
            i = obj.iteration + 1;
            obj.iteration = i;
            
            % Update new state
            x_prev = obj.states(:,i);
            u = obj.controlTraj(:,i);
            x = predictState(obj,x_prev,u);
            obj.states(:,i+1) = x;
            
            % Add any new map features from measurements
            measuredFeatures = obj.measCorrespond(:,i);
            [~,im,~] = intersect(measuredFeatures,obj.mapPointIdx,'stable');
            measuredFeatures(im) = [];
            
            % Add new measured features
            if ~isempty(measuredFeatures)
                obj.mapPointIdx = [obj.mapPointIdx,measuredFeatures];
                obj.map = [obj.map,rand(2,length(measuredFeatures))];
            end
            
            
            % Starting with 1st iteration, re-linearize, re-factor, and
            % re-order sqrt info matrix
            if rem(obj.iteration,100) == 1
                [A,b] = buildSqrtInfoMatrix(obj);
                if nnz(A) > 0
                    p = colamd(A); % re-order A before factoring
                    % update indices
                    if ~isempty(obj.anchorIdx)
                        oldIdx = obj.anchorIdx(:,obj.number);
                        obj.anchorIdx(:,obj.number) = reshape(p(oldIdx),3,1); 
                    end
                    if ~isempty(obj.stateIdx)
                        oldIdx = obj.stateIdx;
                        obj.stateIdx = reshape(p(oldIdx),3,size(oldIdx,2)); 
                    end
                    if ~isempty(obj.mapIdx)
                        oldIdx = obj.mapIdx;
                        obj.mapIdx = reshape(p(oldIdx),2,size(obj.mapIdx,2));
                    end
                    
                    % Re-order A and re-factor
                    A = A(:,p);
%                     b = b(p);
                    [obj.sqrtInfoR, obj.rhsD] = samFactor(A,b);
                end
            else
                % Get rows to add
                [R_add,d_add] = getRowsForSqrtInfoMatrix(obj);
                [obj.sqrtInfoR,obj.rhsD] = incrementalQR(obj.sqrtInfoR,obj.rhsD,R_add,d_add); % update factorization with added rows
            end
            
            % Solve for new estimates
            offset = obj.sqrtInfoR \ obj.rhsD;
            if ~isempty(obj.anchorIdx)
                obj.anchor(:,obj.number) = obj.anchor(:,obj.number) +...
                        reshape(offset(obj.anchorIdx(:,obj.number)),3,1); 
            end
            if ~isempty(obj.stateIdx)
                numStates = size(obj.states,2);
                obj.states = obj.states +...
                        reshape(offset(obj.stateIdx(:)),3,numStates); 
            end
            if ~isempty(obj.mapIdx)
                numMap = size(obj.map,2);
                obj.map = obj.map +...
                        reshape(offset(obj.mapIdx(:)),2,numMap); 
            end
            
        end
        
        % function to plot robot position
        function fig = plotRobot(obj,fig)
%             pos = horzcat(obj.position);
            if nargin < 2
                fig = figure;
            end
            
            % Add Anchor pose for plotting
            statesToPlot = addPose(obj.anchor(:,obj.number),obj.states);
            mapToPlot = addPose(obj.anchor(:,obj.number),obj.map);
            
            hold on;
            plotPose(statesToPlot(1,end),statesToPlot(2,end),statesToPlot(3,end),fig)
            plot(statesToPlot(1,:),statesToPlot(2,:),'b.-');
            plot(mapToPlot(1,:),mapToPlot(2,:),'.');
%             hold off;
        end
        
        
        % function to build full sqrt info matrix
        function [A,b] = buildSqrtInfoMatrix(obj)
            numTimeSteps = size(obj.states,2)-1;
            numTrans = 3*numTimeSteps;
            numStates = numTrans+3;
            measPerTime = size(obj.measurements,1);
            numMeas = measPerTime*numTimeSteps;
            numRows = numTrans + numMeas;
            numMap  = 2*size(obj.map,2);
            numCols = numStates + numMap;
            if obj.hasPrior
                numRows = numRows + 3;
                numCols = numCols + 3;
            end
            A = spalloc(numRows,numCols,1);
            b = zeros(numRows,1);
            
            % Count added variables for indexing
            varCount = 0;
            measureCount = 0;
            G = eye(3); % identity jacobian
            
            if obj.hasPrior
                obj.anchorIdx = [1:3]'; % affected cols
                affectedRows = (measureCount+1):(measureCount+3);
                A(affectedRows,obj.anchorIdx) = obj.invSqrtPriorCov*G;
                varCount = varCount + 3;
                measureCount = measureCount + 3;
                
                % Calculate measurement
                b(affectedRows) = obj.invSqrtPriorCov*(obj.anchor - obj.priorMean);
            end
            
            % Add first state
            obj.stateIdx = zeros(3,numTimeSteps+1);
            newMapEntries = numMap/2 - size(obj.mapIdx,2);
            obj.mapIdx = [obj.mapIdx, zeros(2,newMapEntries)];
            
            obj.stateIdx(:,1) = [(varCount+1):(varCount+3)]';
            varCount = varCount + 3;
            
            % Add other states
            for i = 1:numTimeSteps
                affectedRows = (measureCount+1):(measureCount+3);
                affectedCols = (varCount-2):(varCount+3); % jacobians for prev and current state
                varCount = varCount + 3;
                measureCount = measureCount + 3;
                x_prev = obj.states(:,i);
                x = obj.states(:,i+1);
                u = obj.controlTraj(:,i);                
                
                xp = predictState(obj,x,u);        
                F = getProcessJacobian(obj,x_prev,u);
                
                A(affectedRows,affectedCols) = obj.invSqrtQ*[F,G];
                b(affectedRows) = obj.invSqrtQ*(x - xp);
                
                obj.stateIdx(:,i+1) = affectedCols(4:6)';
            end
            
            % Add measurements
            for i = 1:numTimeSteps
                for j = 1:(measPerTime/2)
                    affectedRows = (measureCount+1):(measureCount+2);
                    measureCount = measureCount + 2;
                    
                    % Two measurements per map point
                    rangeRow = 2*j-1;
                    bearingRow = 2*j;
                    mapPointNum = obj.measCorrespond(j,i);
                    mapCol = find(obj.mapPointIdx == mapPointNum);
                    
                    % Add map location as new var if not added yet                    
                    if any(obj.mapIdx(:,mapCol) == 0)
                        affectedCols = [(varCount+1):(varCount+2)]'; 
                        varCount = varCount + 2;
                        obj.mapIdx(:,mapCol) = affectedCols;
                    else
                        affectedCols = obj.mapIdx(:,mapCol);
                    end
                    % Add state to affect cols
                    affectedCols = [obj.stateIdx(:,i+1);affectedCols];
                    
                    % Get measurement jacobians
                    x = obj.states(:,i+1);
                    m = obj.map(:,mapCol);
                    y = obj.measurements([rangeRow,bearingRow],i);
                    yp = predictMeasurement(obj,x,m);
                    [H,J] = getMeasurementJacobians(obj,x,m);
                    
                    % store in A and b
                    A(affectedRows,affectedCols) = obj.invSqrtR*[H,J];
                    b(affectedRows) = obj.invSqrtR*(y - yp);
                end
            end
            
        end
        
        % Function to add to sqrt info matrix before incremental factoring
        function [R_add,d_add] = getRowsForSqrtInfoMatrix(obj)
            i = obj.iteration;
            numTimeSteps = size(obj.states,2)-1;
            numTrans = 3*numTimeSteps;
            numStates = numTrans+3;
            measPerTime = size(obj.measurements,1);
            numMap = 2*size(obj.map,2);
            numCols = numStates + numMap;
            if obj.hasPrior
                numCols = numCols + 3;
            end
            
            % Pre-allocate
            numRows = 3 + 2*measPerTime;
            R_add = sparse(numRows,numCols);
            d_add = zeros(numRows,1);
            
            measureCount = 0; % don't add rows here, so start at zero     
            varCount = size(obj.sqrtInfoR,2);
            G = eye(3); % identity jacobian
            
            % Add transition constraint            
            affectedRows = (measureCount+1):(measureCount+3);
            affectedCols = [obj.stateIdx(:,i);[(varCount+1):(varCount+3)]']; % jacobians for prev and current state
            varCount = varCount + 3;
            measureCount = measureCount + 3;
            x_prev = obj.states(:,i);
            x = obj.states(:,i+1);
            u = obj.controlTraj(:,i);                

            xp = predictState(obj,x,u);        
            F = getProcessJacobian(obj,x_prev,u);

            R_add(affectedRows,affectedCols) = obj.invSqrtQ*[F,G];
            d_add(affectedRows) = obj.invSqrtQ*(x - xp);
            
            obj.stateIdx(:,i+1) = affectedCols(4:6)';
            
            % Add measurement constraints
            for j = 1:(measPerTime/2)
                affectedRows = (measureCount+1):(measureCount+2);
                measureCount = measureCount + 2;
                    
                % Two measurements per map point
                rangeRow = 2*j-1;
                bearingRow = 2*j;
                mapPointNum = obj.measCorrespond(i,j);
                mapCol = find(obj.mapPointIdx == mapPointNum);

                % Add map location as new var if not added yet                    
                if any(obj.mapIdx(:,mapCol) == 0)
                    affectedCols = (varCount+1):(varCount+2); 
                    varCount = varCount + 2;
                    obj.mapIdx(:,mapCol) = affectedCols';
                else
                    affectedCols = obj.mapIdx(:,mapCol);
                end
                % Add state to affect cols
                affectedCols = [obj.stateIdx(:,i+1);affectedCols];

                % Get measurement jacobians
                x = obj.states(:,i+1);
                m = obj.map(:,mapCol);
                y = obj.measurements([rangeRow,bearingRow],i);
                yp = predictMeasurement(obj,x,m);
                [H,J] = getMeasurementJacobians(obj,x,m);

                % store in A and b
                R_add(affectedRows,affectedCols) = obj.invSqrtR*[H,J];
                d_add(affectedRows) = obj.invSqrtR*(y - yp);
            end
        end
            
        
    end
end