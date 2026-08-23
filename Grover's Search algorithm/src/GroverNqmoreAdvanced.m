clc;
clear;
close all;

%% ==== Get user input ====
n = input('Enter number of qubits: ');

valid = false;
while ~valid
    target_bin = input('Enter target state as a binary string: ', 's');
    if length(target_bin) == n && all(ismember(target_bin, '01'))
        valid = true;
    else
        fprintf('Invalid input. Must be a binary string of length %d.\n', n);
    end
end

N = 2^n;
target = bin2dec(target_bin) + 1;   % convert binary string to index (1-based)
fprintf('\nTarget entered: %s\n', target_bin);
fprintf('Decimal value: %d\n', bin2dec(target_bin));
fprintf('MATLAB index: %d\n', target);
fprintf('Corresponding state: %s\n', dec2bin(target-1,n));
%% ==== Build generalized Hadamard H^{⊗n} ====
H = (1/sqrt(2))*[1 1; 1 -1];

Hn = 1;
for k = 1:n
    Hn = kron(Hn, H);
end

%% ==== Initial state |00...0> ====
psi0 = zeros(N,1);
psi0(1) = 1;

%% ==== Uniform superposition ====
psi = Hn * psi0;

%% ==== Build Oracle ====
w = zeros(N,1);
w(target) = 1;
O = eye(N) - 2*(w*w');

%% ==== Build Diffusion Operator ====
s = ones(N,1)/sqrt(N);
D = 2*(s*s') - eye(N);

%% ==== Optimal number of Grover iterations ====
r = round((pi/4)*sqrt(N));
fprintf('Number of Grover iterations: %d\n', r);
fprintf('Oracle marked state = %s\n', ...
dec2bin(find(diag(O)==-1)-1,n));
%% ==== Run Grover's algorithm ====
for iter = 1:r
    psi = O * psi;      % apply Oracle
    psi = D * psi;       % apply Diffusion
end

%% ==== Compute probabilities ====
P = abs(psi).^2;

%% ==== Display results ====
labels = cell(1,N);
for k = 0:N-1
    labels{k+1} = dec2bin(k, n);
end
x=1:N;
figure;
bar(x,P);
set(gca,'XTick',x, 'XTickLabel', labels);
xlabel('Basis State');
ylabel('Probability');
title(['Grover Search — n = ', num2str(n), ', Target = |', target_bin, '>']);

fprintf('\nFinal probabilities:\n')
for k = 1:N
    fprintf('%s : %.4f\n', labels{k}, P(k))
end