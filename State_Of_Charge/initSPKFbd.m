% init SPKF for bar-delta method
function spkfData = initSPKFbd(v0, T0, SigmaX0, SigmaV, SigmaW, model)
  % First, initialize the variables for the "bar" filter
  % Initial state description
  ir0 = 0; spkfData.irInd = 1;
  hk0 = 0; spkfData.hkInd = 2;
  SOC0 = mean(SOCfromOCVtemp(v0, T0, model)); spkfData.zkInd = 3;
  ib0 = 0; spkfData.biasInd = 4; % Add state to estimate current sensor bias
  spkfData.xhat = [ ir0 hk0 SOC0 ib0]'; % initial state

  % Covariance values
  spkfData.SigmaX = SigmaX0;
  spkfData.SigmaV = SigmaV;
  spkfData.SigmaW = SigmaW;
  spkfData.Snoise = real(chol(blkdiag(SigmaW, SigmaV), 'lower' ));
  spkfData.Qbump = 5;

  % SPKF specific parameters
  Nx = length(spkfData.xhat); spkfData.Nx = Nx; % state-vector length
  Ny = 1; spkfData.Ny = Ny; % measurement-vector length
  Nu = 1; spkfData.Nu = Nu; % input-vector length
  Nw = size(SigmaW,1); spkfData.Nw = Nw; % process-noise-vector length
  Nv = size(SigmaV,1); spkfData.Nv = Nv; % sensor-noise-vector length
  Na = Nx+Nw+Nv; spkfData.Na = Na; % augmented-state-vector length

  h = sqrt(3); spkfData.h = h; % SPKF/CDKF tuning factor
  Weight1 = (h*h-Na)/(h*h); % weighting factors when computing mean
  Weight2 = 1/(2*h*h); % and covariance
  spkfData.Wm = [ Weight1; Weight2*ones(2*Na,1)]; % mean
  spkfData.Wc = spkfData.Wm; % covar

  % previous value of current
  spkfData.priorI = 0;
  spkfData.signIk = 0;

  % store model data structure too
  spkfData.model = model;
  
  % Now, initialize variables for the "delta" filters
  % SOC is estimated using SPKF
  dNx = 1; spkfData.dNx = dNx; % one state per delta filter, for estimating SOC
  dNw = 1; spkfData.dNw = dNw; % one process-noise per delta filter
  dNa = dNx+dNw; spkfData.dNa = dNa; % augmented length for delta filters
  
  spkfData.dh = h; 
  Weight1 = (h*h-dNa)/(h*h); % weighting factors when computing mean
  Weight2 = 1/(2*h*h); % and covariance
  spkfData.dWm = [ Weight1; Weight2*ones(2*dNa,1)]; % mean
  spkfData.dWc = spkfData.dWm; % covar
  
  spkfData.celldz = 0*v0(:);
  spkfData.cellSdz = SigmaX0(spkfData.zkInd,spkfData.zkInd)*ones(size(v0(:)));  
end