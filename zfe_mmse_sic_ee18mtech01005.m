          %%% PRABHAT KUMAR RAI %%%
          %%%  EE18MTECH01005  %%%
%%% Performance of Spatial Multiplexing %%%

clear all; close all; clc;
itr = 100000;
N = 12; % number of bits or symbols
Eb_N0_dB = 12:3:24; % multiple SNR values
% Eb_N0_dB = 19.5;
nt = 4; nr = 6;

for p = 1:length(Eb_N0_dB)
    num_error_zfe_sic = 0; num_error_mmse = 0; num_error_zfe = 0; num_error_mmse_sic = 0;
    
    for ii = 1 : itr
        input_1 = rand(1, N/nt)>0.5; % generating 0,1 with equal probability
        input_2 = rand(1, N/nt)>0.5; % generating 0,1 with equal probability
        input_3 = rand(1, N/nt)>0.5; % generating 0,1 with equal probability
        input_4 = rand(1, N/nt)>0.5; % generating 0,1 with equal probability

        %%% PSK gray coding
        psk_comp_0 = 1+0i; psk_comp_4 = -1+0i;
        psk_comp_1 = (1+1i)/sqrt(2); psk_comp_5 = (-1-1i)/sqrt(2);
        psk_comp_2 = 0+1i; psk_comp_6 = 0-1i;
        psk_comp_3 = (-1+1i)/sqrt(2); psk_comp_7 = (1-1i)/sqrt(2);  
        
        x1 = mpsk(input_1); x2 = mpsk(input_2);
        x3 = mpsk(input_3); x4 = mpsk(input_4);
        input = [input_1 input_2 input_3 input_4];
        
        %%% 8-PSK signal to transmit  
        signal_psk = [x1; x2; x3; x4];
        signal_psk = signal_psk;

        %%% Channel
        channel = (rand(nr, nt)+1i*rand(nr, nt));
        signal = channel*signal_psk;
        
        %%% Signal Reception with Noise Addition
        noise = (1/sqrt(2))*(randn(nr, 1) + 1i*randn(nr, 1));
        receive_signal = signal + 10^(-Eb_N0_dB(p)/20)*noise;
        
        
%%%%%%%%%%%%--------- ZFE-SIC equalization based Detection -------%%%%%%
        
        [Q, R] = qr(channel, 0);    %% Q-R decomposition
        
        y_tilda = Q'*receive_signal; 

        x_hat_final_zfe_sic = []; output_zfe_sic = []; x_hat4_zfe_sic = 0; x_hat3_zfe_sic = 0; x_hat2_zfe_sic = 0;
        
        for ll = 0 : (nt-1)
            minD_0 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_0);
            minD_1 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_1);
            minD_2 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_2);
            minD_3 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_3);
            minD_4 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_4);
            minD_5 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_5);
            minD_6 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_6);
            minD_7 = abs(y_tilda(nt-ll) - x_hat4_zfe_sic*R(nt-ll, nt) - x_hat3_zfe_sic*R(nt-ll, nt-1) - x_hat2_zfe_sic*R(nt-ll, nt-2) - R(nt-ll, nt-ll)*psk_comp_7);
            
            if (minD_0 < minD_1 && minD_0 < minD_2 && minD_0 < minD_3 && minD_0 < minD_4 && minD_0 < minD_5 && minD_0 < minD_6 && minD_0 < minD_7)
                output_zfe_sic = [ 0 0 0 output_zfe_sic ]; x_hat_final_zfe_sic = [ 1+0i x_hat_final_zfe_sic ];
            elseif (minD_1 < minD_0 && minD_1 < minD_2 && minD_1 < minD_3 && minD_1 < minD_4 && minD_1 < minD_5 && minD_1 < minD_6 && minD_1 < minD_7)
                output_zfe_sic = [ 0 0 1 output_zfe_sic ]; x_hat_final_zfe_sic = [ (1+1i)/sqrt(2) x_hat_final_zfe_sic ];   
            elseif (minD_2 < minD_0 && minD_2 < minD_1 && minD_2 < minD_3 && minD_2 < minD_4 && minD_2 < minD_5 && minD_2 < minD_6 && minD_2 < minD_7)
                output_zfe_sic = [ 1 0 1 output_zfe_sic ]; x_hat_final_zfe_sic = [ 0+1i x_hat_final_zfe_sic ];
            elseif (minD_3 < minD_0 && minD_3 < minD_1 && minD_3 < minD_2 && minD_3 < minD_4 && minD_3 < minD_5 && minD_3 < minD_6 && minD_3 < minD_7)
                output_zfe_sic = [ 1 0 0 output_zfe_sic ]; x_hat_final_zfe_sic = [ (-1+1i)/sqrt(2) x_hat_final_zfe_sic ];                
            elseif (minD_4 < minD_0 && minD_4 < minD_2 && minD_4 < minD_3 && minD_4 < minD_1 && minD_4 < minD_5 && minD_4 < minD_6 && minD_4 < minD_7)            
                output_zfe_sic = [ 1 1 0 output_zfe_sic ]; x_hat_final_zfe_sic = [ -1+0i x_hat_final_zfe_sic ];                
            elseif (minD_5 < minD_0 && minD_5 < minD_2 && minD_5 < minD_3 && minD_5 < minD_4 && minD_5 < minD_1 && minD_5 < minD_6 && minD_5 < minD_7)            
                output_zfe_sic = [ 1 1 1 output_zfe_sic ]; x_hat_final_zfe_sic = [ (-1-1i)/sqrt(2) x_hat_final_zfe_sic ];                
            elseif (minD_6 < minD_0 && minD_6 < minD_2 && minD_6 < minD_3 && minD_6 < minD_4 && minD_6 < minD_5 && minD_6 < minD_1 && minD_6 < minD_7)            
                output_zfe_sic = [ 0 1 1 output_zfe_sic ]; x_hat_final_zfe_sic = [ 0-1i x_hat_final_zfe_sic ];                  
            elseif (minD_7 < minD_0 && minD_7 < minD_2 && minD_7 < minD_3 && minD_7 < minD_4 && minD_7 < minD_5 && minD_7 < minD_6 && minD_7 < minD_1)
                output_zfe_sic = [ 0 1 0 output_zfe_sic ]; x_hat_final_zfe_sic = [ (1-1i)/sqrt(2) x_hat_final_zfe_sic ]; 
            end
            x_hat4_zfe_sic = x_hat_final_zfe_sic(ll+1);
            if (ll>=1)
                x_hat3_zfe_sic = x_hat_final_zfe_sic(ll);
            else
                x_hat3_zfe_sic = 0;
            end
            
            if (ll>=2)
                x_hat2_zfe_sic = x_hat_final_zfe_sic(ll-1);
            else
                x_hat2_zfe_sic = 0;
            end
        end
        
        % Error for ZFE-SIC
        num_error_zfe_sic = num_error_zfe_sic + (size(find([input - output_zfe_sic]),2)); 
        
%%%%%%%%%%%%%%--------- MMSE equalization based Detection -------%%%%%%
       
        SNR = 10^(Eb_N0_dB(p)+10*log10(3));
        
        B_mmse = ((channel'*channel) + (nt/SNR)*eye(nt))\(channel');
        
        recev_mmse = B_mmse*receive_signal;
        
        x_hat_final_mmse = []; output_mmse = [];
        
        for ll = 1 : nt
            minD_0 = norm(recev_mmse(ll) - psk_comp_0);
            minD_1 = norm(recev_mmse(ll) - psk_comp_1);
            minD_2 = norm(recev_mmse(ll) - psk_comp_2);
            minD_3 = norm(recev_mmse(ll) - psk_comp_3);
            minD_4 = norm(recev_mmse(ll) - psk_comp_4);
            minD_5 = norm(recev_mmse(ll) - psk_comp_5);
            minD_6 = norm(recev_mmse(ll) - psk_comp_6);
            minD_7 = norm(recev_mmse(ll) - psk_comp_7);
            
            if (minD_0 < minD_1 && minD_0 < minD_2 && minD_0 < minD_3 && minD_0 < minD_4 && minD_0 < minD_5 && minD_0 < minD_6 && minD_0 < minD_7)
                output_mmse = [ output_mmse 0 0 0 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_0 ];
            elseif (minD_1 < minD_0 && minD_1 < minD_2 && minD_1 < minD_3 && minD_1 < minD_4 && minD_1 < minD_5 && minD_1 < minD_6 && minD_1 < minD_7)
                output_mmse = [ output_mmse 0 0 1 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_1 ];   
            elseif (minD_2 < minD_0 && minD_2 < minD_1 && minD_2 < minD_3 && minD_2 < minD_4 && minD_2 < minD_5 && minD_2 < minD_6 && minD_2 < minD_7)
                output_mmse = [ output_mmse 1 0 1 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_2 ];
            elseif (minD_3 < minD_0 && minD_3 < minD_1 && minD_3 < minD_2 && minD_3 < minD_4 && minD_3 < minD_5 && minD_3 < minD_6 && minD_3 < minD_7)
                output_mmse = [ output_mmse 1 0 0 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_3 ];                
            elseif (minD_4 < minD_0 && minD_4 < minD_2 && minD_4 < minD_3 && minD_4 < minD_1 && minD_4 < minD_5 && minD_4 < minD_6 && minD_4 < minD_7)            
                output_mmse = [ output_mmse 1 1 0 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_4 ];                
            elseif (minD_5 < minD_0 && minD_5 < minD_2 && minD_5 < minD_3 && minD_5 < minD_4 && minD_5 < minD_1 && minD_5 < minD_6 && minD_5 < minD_7)            
                output_mmse = [ output_mmse 1 1 1 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_5 ];                
            elseif (minD_6 < minD_0 && minD_6 < minD_2 && minD_6 < minD_3 && minD_6 < minD_4 && minD_6 < minD_5 && minD_6 < minD_1 && minD_6 < minD_7)            
                output_mmse = [ output_mmse 0 1 1 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_6 ];                  
            elseif (minD_7 < minD_0 && minD_7 < minD_2 && minD_7 < minD_3 && minD_7 < minD_4 && minD_7 < minD_5 && minD_7 < minD_6 && minD_7 < minD_1)
                output_mmse = [ output_mmse 0 1 0 ]; x_hat_final_mmse = [ x_hat_final_mmse psk_comp_7 ]; 
            end
        end
        % Error for MMSE
        num_error_mmse = num_error_mmse + (size(find([input - output_mmse]),2)); 
        
%%%%%%%%%%%%%%--------- MMSE-SIC equalization based Detection -------%%%%%%
    x_hat_mmse_sic = []; output_mmse_sic = []; x_hat_final_mmse_sic=[];
    
    for ll = 0 : nt-1
        SNR = 10^(Eb_N0_dB(p)+10*log10(3));
        B_mmse_sic = ((channel(:,1:nt-ll)'*channel(:,1:nt-ll)) + ((nt-ll)/SNR)*eye(nt-ll))\(channel(:,1:nt-ll)');
        recev_mmse_sic_1 = B_mmse_sic(nt-ll:nt-ll:end)*receive_signal;
        
        minD_0 = norm(recev_mmse_sic_1 - psk_comp_0);
        minD_1 = norm(recev_mmse_sic_1 - psk_comp_1);
        minD_2 = norm(recev_mmse_sic_1 - psk_comp_2);
        minD_3 = norm(recev_mmse_sic_1 - psk_comp_3);
        minD_4 = norm(recev_mmse_sic_1 - psk_comp_4);
        minD_5 = norm(recev_mmse_sic_1 - psk_comp_5);
        minD_6 = norm(recev_mmse_sic_1 - psk_comp_6);
        minD_7 = norm(recev_mmse_sic_1 - psk_comp_7);
        
        if (minD_0 < minD_1 && minD_0 < minD_2 && minD_0 < minD_3 && minD_0 < minD_4 && minD_0 < minD_5 && minD_0 < minD_6 && minD_0 < minD_7)
            output_mmse_sic = [ 0 0 0 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_0;
        elseif (minD_1 < minD_0 && minD_1 < minD_2 && minD_1 < minD_3 && minD_1 < minD_4 && minD_1 < minD_5 && minD_1 < minD_6 && minD_1 < minD_7)
            output_mmse_sic = [ 0 0 1 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_1;       
        elseif (minD_2 < minD_0 && minD_2 < minD_1 && minD_2 < minD_3 && minD_2 < minD_4 && minD_2 < minD_5 && minD_2 < minD_6 && minD_2 < minD_7)
            output_mmse_sic = [ 1 0 1 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_2;  
        elseif (minD_3 < minD_0 && minD_3 < minD_1 && minD_3 < minD_2 && minD_3 < minD_4 && minD_3 < minD_5 && minD_3 < minD_6 && minD_3 < minD_7)
            output_mmse_sic = [ 1 0 0 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_3;                
        elseif (minD_4 < minD_0 && minD_4 < minD_2 && minD_4 < minD_3 && minD_4 < minD_1 && minD_4 < minD_5 && minD_4 < minD_6 && minD_4 < minD_7)            
            output_mmse_sic = [ 1 1 0 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_4;                 
        elseif (minD_5 < minD_0 && minD_5 < minD_2 && minD_5 < minD_3 && minD_5 < minD_4 && minD_5 < minD_1 && minD_5 < minD_6 && minD_5 < minD_7)            
            output_mmse_sic = [ 1 1 1 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_5;                
        elseif (minD_6 < minD_0 && minD_6 < minD_2 && minD_6 < minD_3 && minD_6 < minD_4 && minD_6 < minD_5 && minD_6 < minD_1 && minD_6 < minD_7)            
            output_mmse_sic = [ 0 1 1 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_6;                    
        elseif (minD_7 < minD_0 && minD_7 < minD_2 && minD_7 < minD_3 && minD_7 < minD_4 && minD_7 < minD_5 && minD_7 < minD_6 && minD_7 < minD_1)
            output_mmse_sic = [ 0 1 0 output_mmse_sic ]; x_hat_final_mmse_sic = psk_comp_7; 
        end
        x_hat_mmse_sic = [x_hat_final_mmse_sic x_hat_mmse_sic ];
        receive_signal = receive_signal - x_hat_final_mmse_sic*channel(:,nt-ll);      
    end
    
    % Error for MMSE-SIC
    num_error_mmse_sic = num_error_mmse_sic + (size(find([input - output_mmse_sic]),2)); 
        
    end
    ber_zfe_sic(p) = num_error_zfe_sic/(N*itr); % simulated ber
    ber_mmse(p) = num_error_mmse/(N*itr); % simulated ber
    ber_mmse_sic(p) = num_error_mmse_sic/(N*itr); % simulated ber
end
% ber_zfe_sic(end) = 7.1*10^-4;
% ber_mmse(end) = 1.25*10^-3;
% ber_mmse_sic(end) = 6.3*10^-4;
figure;
semilogy(Eb_N0_dB, ber_zfe_sic,'mo-','LineWidth',2); 
hold on; 
semilogy(Eb_N0_dB, ber_mmse,'b-','LineWidth',2); 
semilogy(Eb_N0_dB, ber_mmse_sic,'g-','LineWidth',2); 
axis([Eb_N0_dB(1) Eb_N0_dB(end) 10^-5 0.5]); grid on;
legend('ZFE-SIC','MMSE','MMSE-SIC'); xlabel('SNR in dB'); ylabel('Probability of error, Pe');
title('Performance of Spatial Multiplexing by EE18MTECH01005');


function minD = mpsk(x)
         psk_0 = [0 0 0]; psk_4 = [1 1 0];
         psk_1 = [0 0 1]; psk_5 = [1 1 1];
         psk_2 = [1 0 1]; psk_6 = [0 1 1];
         psk_3 = [1 0 0]; psk_7 = [0 1 0];
         if (x == psk_0)
             minD = 1+0i;
         elseif (x == psk_1)
             minD = (1+1i)/sqrt(2);  
         elseif (x == psk_2)
             minD = 0+1i;    
         elseif (x == psk_3)
             minD = (-1+1i)/sqrt(2);
         elseif (x == psk_4)
             minD = -1+0i;    
         elseif (x == psk_5)
             minD = (-1-1i)/sqrt(2);    
         elseif (x == psk_6)
             minD = 0-1i;    
         elseif (x == psk_7)
             minD = (1-1i)/sqrt(2);        
         end
end