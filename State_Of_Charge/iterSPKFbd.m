function [zkbar, zkbarbnd, dzk, dzkbnd, ib, ibbnd, spkfData] = ...
         iterSPKFbd(vk, ik, Tk, deltat, spkfData)
  model = spkfData. model;

  % Load the cell model parameters
  G = getParamESC('GParam' , Tk, model);
  M = getParamESC('MParam' , Tk, model);
  RC = exp(-deltat./abs(getParamESC('RCParam' , Tk, model)))';
  R = getParamESC('RParam' , Tk, model)';
  R0 = getParamESC('R0Param' , Tk, model);

  % Get data stored in spkfData structure
  I = spkfData.priorI;
  SigmaX = spkfData.SigmaX;
  xhat = spkfData.xhat;
  celldz = spkfData.celldz;
  cellSdz = spkfData.cellSdz;
  Nx = spkfData.Nx;
  Nw = spkfData.Nw;
  Nv = spkfData.Nv;
  Na = spkfData.Na;
  Snoise = spkfData.Snoise;
  Wc = spkfData.Wc;
  irInd = spkfData.irInd;
  hkInd = spkfData.hkInd;
  zkInd = spkfData.zkInd;
  ibInd = spkfData.biasInd; % Add index for current-sensor bias estimator
  
  dQinv = spkfData.dQinv;
  Qinvbar = spkfData.Qinvbar;

  % Step 1a : State estimate time update
  %       - Create xhatminus augmented SigmaX points
  %       - Extract xhatminus state SigmaX points
  %       - Compute weighted average xhatminus(k)

  % Step 1a-1: Create augmented SigmaX and xhat
  [ sigmaXa, p] = chol(SigmaX, 'lower' );
  if p>0,
    fprintf('Cholesky error. Recovering...\n' );
    theAbsDiag = abs(diag(SigmaX));
    sigmaXa = diag(max(SQRT(theAbsDiag), SQRT(spkfData.SigmaW(1,1))));
  end
  sigmaXa=[real(sigmaXa) zeros([Nx Nw+Nv]); zeros([Nw+Nv Nx]) Snoise];
  xhata = [ xhat; zeros([Nw+Nv 1])];
  % NOTE: sigmaX a is lower-triangular

  % Step 1a-2: Calculate SigmaX points (strange indexing of xhat a to
  % avoid "repmat" call, which is very inefficient in MATLAB)
  Xa = xhata(:, ones([1 2*Na+1])) + ...
      spkfData.h*[zeros([ Na 1]), sigmaXa, -sigmaXa];

  % Step 1a-3: Time update from last iteration until now
  %   state Eqn(xold, current, xnoise)
  Xx = stateEqn(Xa(1: Nx,:), I, Xa(Nx+1: Nx+Nw,:));
  xhat = Xx*spkfData.Wm;

  % Step 1b: Error covariance time update
  %           - Compute weighted covariance sigmaminus(k)
  %           (strange indexing of xhat to avoid "repmat" call)
  Xs = Xx - xhat(:,ones([1 2*Na+1]));
  SigmaX = Xs*diag(Wc)*Xs';

  % Step 1c: Output estimate
  %           - Compute weighted output estimate yhat(k)
  I = ik; yk = vk;
  Y = outputEqn(Xx, I, Xa(Nx+Nw+1: end,:), Tk, model);
  yhat = Y*spkfData.Wm;

  % Step 2a: Estimator gain matrix
  Ys = Y - yhat(:, ones([1 2*Na+1]));
  SigmaXY = Xs*diag(Wc)*Ys';
  SigmaY = Ys*diag(Wc)*Ys';
  L = SigmaXY/SigmaY;

  % Step 2b: State estimate measurement update
  r = mean(yk) - yhat; % residual. Use to check for sensor errors...
  if r^2 > 100*SigmaY, L(:,1)=0.0; end
  xhat = xhat + L*r;
  xhat(zkInd)=min(1.05, max(-0.05, xhat(zkInd)));

  % Step 2c: Error covariance measurement update
  SigmaX = SigmaX - L*SigmaY*L';
  [~, S, V] = svd(SigmaX);
  HH = V*S*V';
  SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness

  % Q-bump code
  if r^2>4*SigmaY, % bad voltage estimate by 2-SigmaX, bump Q
      fprintf('Bumping sigmax\n' );
      SigmaX(zkInd, zkInd) = SigmaX(zkInd, zkInd)*spkfData.Qbump;
  end

  % Save data in spkfData structure for next time...
  ib = xhat(ibInd);
  ibbnd = 3*sqrt(SigmaX(ibInd, ibInd));
  spkfData.priorI = ik;
  spkfData.SigmaX = SigmaX;
  spkfData.xhat = xhat;
  zkbar = xhat(zkInd);
  zkbarbnd = 3*sqrt(SigmaX(zkInd, zkInd));
  
  % The "bar" filter update is complete. Now, work on "delta" filter updates
  offset = M*xhat(hkInd) - R*xhat(irInd) - R0*(ik-ib);
  for thecell = 1:length(vk),
    % Implement SPKF for delta-soc
    I = spkfData.priorI;
    % First, update the delta-SOC SPKFs
    % Step 1a - State prediction time update
    cellxa = [celldz(thecell); 0];
    cellSxa = diag([SQRT(cellSdz(thecell)), sqrt(spkfData.SigmaW(1,1))]);
    cellXa = cellxa(:,ones([1 2*spkfData.dNa+1])) + ...
             spkfData.h*[0*cellxa,cellSxa,-cellSxa];
    cellXx = cellStateEqn(cellXa,I-ib,dQinv(thecell));
    celldz(thecell) = cellXx*spkfData.dWm;
    % Step 1b - Do error covariance time update
    cellXs = cellXx - celldz(thecell);
    cellSdz(thecell) = cellXs*diag(spkfData.dWc)*cellXs';
    % Step 1c - output estimate
    I = ik;     
    cellY = cellOutputEqn(cellXx,I-ib,cellXa,zkbar,spkfData.dR0(thecell),offset,Tk,model);
    cellyhat = cellY*spkfData.dWm;
    % Step 2a - Estimator gain matrix
    cellYs = cellY - cellyhat;
    cellSy = cellYs*diag(spkfData.dWc)*cellYs' + spkfData.SigmaV; % lin sensor noise
    cellSxy = cellXs*diag(spkfData.dWc)*cellYs';
    cellL = cellSxy/cellSy;
    % Step 2b - State estimate measurement update
    celldz(thecell) = celldz(thecell) + cellL*(vk(thecell) - cellyhat);
    % Step 2c - Error covariance measurement update
    cellSdz(thecell) = cellSdz(thecell) - cellL*cellSy*cellL;
  end

  % Save data in spkfData structure for next time...
  spkfData.celldz = celldz;
  spkfData.cellSdz = cellSdz;
  dzk = celldz;
  dzkbnd = 3*sqrt(cellSdz);

  % Calculate new states for all of the old state vectors in xold.
  function xnew = stateEqn(xold, current, xnoise)
    current = current + xnoise(1,:); % noise adds to current
    xnew = 0*xold;
    xnew(irInd,:) = RC*xold(irInd,:) + (1-RC)*(current-xold(ibInd,:));
    Ah = exp(-abs((current+xold(ibInd,:))*G*deltat*Qinvbar/3600)); % hysteresis factor
    xnew(hkInd,:) = Ah.*xold(hkInd,:) + (Ah-1).*sign(current-xold(ibInd,:));
    xnew(zkInd,:) = xold(zkInd,:) - (current-xold(ibInd,:))*Qinvbar/3600;
    xnew(hkInd,:) = min(1, max(-1, xnew(hkInd,:)));
    xnew(zkInd,:) = min(1.05, max(-0.05, xnew(zkInd,:)));
    xnew(ibInd,:) = xold(ibInd,:) + xnoise(2,:);
  end

  function xnew = cellStateEqn(xold,oldI,celldQinv)
    xnew = xold(1,:) - (oldI + xold(2,:))*celldQinv/3600; % use "known" inverse capacities
  end

  % Calculate cell output voltage for all of state vectors in xhat
  function yhat = outputEqn(xhat, current, ynoise, T, model)
    yhat = OCVfromSOCtemp(xhat(zkInd,:), T, model);
    yhat = yhat + M*xhat(hkInd,:); 
    yhat = yhat - R*xhat(irInd,:) - R0*(current - xhat(ibInd,:)) + ynoise(1,:);
  end

  function yhat = cellOutputEqn(cellXx,I,cellXa,Z,dR,offset,T,model)
    I = I + cellXa(2,:); % add current noise to current
    yhat = OCVfromSOCtemp(Z+cellXx,T,model); % OCV part
    yhat = yhat - I * dR; % delta resistance part
    yhat = yhat + offset; % polarization, hysteresis, bar resistance
    % note: sensor noise handled separately
  end
  
  % "Safe" square root
  function X = SQRT(x)
      X = sqrt(max(0, x));
  end
end