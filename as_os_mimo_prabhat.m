% PRABHAT KUMAR RAI --- EE18MTECH01005
% MIMO ASSIGNMENT-1 ( Ber Curve for Antipodal & Orthogonal Signalling )

clear all; close all; clc;
N = 3000; % number of bits or symbols
tic;
antip_err_total = 0; ortho_err_total = 0;

Itr = 100; % Number of Iterations
for rr = 1 : Itr
    
    input = rand(1, N) > 0.5; % generating 0,1 with equal probability
    
    % Antipodal & Orthogonal Transmitter
    antipodal_input = [;]; antipodal_a = [1;-1]; antipodal_b = [-1;1];  
    orthogonl_input = [;]; orthogonl_a = [1;0]; orthogonl_b = [0;1];

    for pp = 1 : N
        if (input(pp) == 0)
            antipodal_input = [antipodal_input antipodal_a];
            orthogonl_input = [orthogonl_input orthogonl_a];
        else
            antipodal_input = [antipodal_input antipodal_b];
            orthogonl_input = [orthogonl_input orthogonl_b];
        end
    end

    % Addition of noise
    noise1 = (1/sqrt(2))*(randn(1, N)+1i*randn(1, N));
    noise2 = (1/sqrt(2))*(randn(1, N)+1i*randn(1, N));
    noise = [noise1; noise2];

    % or Noise can be generated like this
%     noise = (1/sqrt(2))*(randn(2, N)+1i*randn(2, N));
    
    EnB = 0 : 14;   % multiple Eb/No values in dB's

    for p = 1 : length(EnB)
        
        % Noise addition to Antipodal
        antipodal_corrp = antipodal_input + 10^(-EnB(p)/20)*(sqrt(2))*noise; % AWGN with different dB's
   
        % Noise addition to Orthogonal
        orthogonl_corrp = orthogonl_input + 10^(-EnB(p)/20)*noise; % AWGN with different dB's

        % Demodulating the signal
        antipodal_out = []; orthogonl_out = [];
        for q = 1 : 2 : (2*length(antipodal_corrp))
            
            % Taking first 2 values from a matrix (2*N)
            antipodal_hg = [antipodal_corrp(q); antipodal_corrp(q+1)]; 
            orthogonl_hg = [orthogonl_corrp(q); orthogonl_corrp(q+1)];

            antipodal_dd_1 = norm(antipodal_hg - antipodal_a); % ||Y - Ua||
            orthogonl_dd_1 = norm(orthogonl_hg - orthogonl_a);         % ||Y - Ua||
       
            antipodal_dd_0 = norm(antipodal_hg - antipodal_b); % ||Y - Ub||
            orthogonl_dd_0 = norm(orthogonl_hg - orthogonl_b);         % ||Y - Ub||
     
            if antipodal_dd_0 > antipodal_dd_1   % ||Y-Ub|| >  ||Y-Ua||
                out = 0;                         % Then output = 0
            else                                 % ||Y-Ub|| <  ||Y-Ua||
                out = 1;                         % Then output = 1
            end
            antipodal_out = [antipodal_out out]; 
       
            if orthogonl_dd_0 > orthogonl_dd_1           % ||Y-Ub|| >  ||Y-Ua||
                outt = 0;                        % Then output = 0
            else                                 % ||Y-Ub|| <  ||Y-Ua||
                outt = 1;                        % Then output = 1
            end
            orthogonl_out = [orthogonl_out outt];
        end
            
        % calculating the number of bit errors
        antipodal_nError(p) = size(find([input - antipodal_out]),2);
        orthogonl_nError(p) = size(find([input - orthogonl_out]),2);
    end
      antip_err_total = antip_err_total + antipodal_nError; % Total Error
      ortho_err_total = ortho_err_total + orthogonl_nError; % Total Error
end
antip_err = antip_err_total/(N*Itr); ortho_err = ortho_err_total/(N*Itr);
% theory_Ber = 0.5*erfc(sqrt(10.^(EnB/10))); figure; 
semilogy(EnB, antip_err, 'r*-', 'linewidth', 1.5); hold on; semilogy(EnB, ortho_err, 'bx-', 'linewidth', 1.5); 
% semilogy(EnB, theory_Ber,'b.-');
axis([0 14 10^-6 0.5]); grid on; legend('Antipodal BER', 'Orthogonal BER'); xlabel('Eb/No in dB');
ylabel('Probability of error, Pe'); title('BER curve for Antipodal and Orthogonal Signalling');
toc;
