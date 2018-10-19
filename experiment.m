function reject_null = experiment(str_loop_p, str_loop_c, str_noise)
T = 100; % no. observations for each subject
p = 30; % no. brain regions
N = 10; % no. patients and controls

alpha = 0.01; % for hypothesis testing

var = .5; % variance of each node

n_loop = 5; % no. nodes in loop
%str_loop_p = .1; % strength of connections along loop
%str_loop_c = .1

%str_noise = .1; % strength of connections not in loop

%% Construct inverse covariance matrices
% Loop
invcov_loop_p = eye(p)*var;
for i = 1:n_loop-1
   invcov_loop_p(i, i+1) = str_loop_p;
   invcov_loop_p(i+1, i) = str_loop_p;
end
invcov_loop_p(1, n_loop) = str_loop_p;
invcov_loop_p(n_loop, 1) = str_loop_p;

invcov_loop_c = eye(p)*var;
for i = 1:n_loop-1
   invcov_loop_c(i, i+1) = str_loop_c;
   invcov_loop_c(i+1, i) = str_loop_c;
end
invcov_loop_c(1, n_loop) = str_loop_c;
invcov_loop_c(n_loop, 1) = str_loop_c;

% Add up to 5 additional noise edges per subject
invcov_noise = zeros(p, p, N);
for n = 1:N
    for i = 1:5
        first = ceil(rand * p); second = ceil(rand * p);
        invcov_noise(first, second, n) = str_noise;
        invcov_noise(second, first, n) = str_noise;
    end
    % and erase the ones that overlap with the cycle
    for i = 1:n_loop-1
        invcov_noise(i, i+1, n) = 0;
        invcov_noise(i+1, i, n) = 0;
    end
    invcov_noise(1, n_loop, n) = 0;
    invcov_noise(n_loop, 1, n) = 0;
end

% combine cycle and noise for patients
invcov_patient = invcov_noise;
for n = 1:N
   invcov_patient(:,:,n) = invcov_patient(:,:,n) + invcov_loop_p;
end

% combine cycle and noise for controls
invcov_control = invcov_noise;
for n = 1:N
    invcov_control(:,:,n) = invcov_control(:,:,n) + invcov_loop_c;
end

val_patient = zeros(1,N); val_control = zeros(1,N);
for n = 1:N
    %% Generate simulation data
    % rng(13); % set random number seed
    cov_patient = (inv(invcov_patient(:,:,n)) + inv(invcov_patient(:,:,n)).')/2;
    cov_control = (inv(invcov_control(:,:,n)) + inv(invcov_control(:,:,n)).')/2;
    X_patient = mvnrnd(zeros(1,p), cov_patient, T); % T by p
    X_control = mvnrnd(zeros(1,p), cov_control, T);
    
    %% Hypothesis testing
    T_patient = RunDebiasedGLasso(X_patient, N, p); % p by p
    T_control = RunDebiasedGLasso(X_control, N, p);
    
    T_patient = GetCorrelation(T_patient);
    T_control = GetCorrelation(T_control);
    
    val_patient(n) = GetTestStatistic(T_patient, n_loop);
    val_control(n) = GetTestStatistic(T_control, n_loop);
end

%colormap jet
%imagesc(T_patient)
%colorbar

reject_null = TwoSampleTTest(val_patient, val_control, alpha);
end

function T = RunDebiasedGLasso(X, N, p)
S = (1/N)*(X.'*X);
rho = real(sqrt(log(p)/N));
%[Theta,~] = ADMM(S,rho,1)
[Theta,~] = graphicalLasso(S,rho,100,.01);

% Debiasing
T = 2 * Theta - Theta * S * Theta;

% T = inv(S);
end

function val = GetTestStatistic(T, n_loop)
% add up values around loop in debiased precision matrix T
val = 0;

for i = 1:n_loop-1
   val = val + abs(T(i, i+1));
end
val = val + abs(T(1, n_loop));

end

function coh = GetCorrelation(S_inv)
% normalize off-diagonal entries by diagonal entries
[p, ~] = size(S_inv);
for i=1:p
    for j=1:p
        coh(i,j)=S_inv(i,j)/sqrt(S_inv(i,i)*S_inv(j,j));
    end
end
end