          %%% PRABHAT KUMAR RAI %%%
          %%%  EE18MTECH01005  %%%
%%% Alamouti Code with 16-QAM Constellation %%%


clear all; close all; clc;
itr = 1000000;
N = 8; % number of bits or symbols
Eb_N0_dB = 5:5:35; % multiple SNR values
for p = 1:length(Eb_N0_dB)
    
    bererror = 0; sererror = 0; cererror = 0;
    
    for ii = 1 : itr
        input_1 = rand(1,N/2)>0.5; % generating 0,1 with equal probability
        input_2 = rand(1,N/2)>0.5; % generating 0,1 with equal probability
        input = [input_1 input_2];

        %%% QAM gray coding
        qam_0 = [0 0 0 0]; qam_4 = [0 1 0 0]; qam_8 = [1 1 0 0]; qam_12 = [1 0 0 0];
        qam_1 = [0 0 0 1]; qam_5 = [0 1 0 1]; qam_9 = [1 1 0 1]; qam_13 = [1 0 0 1];
        qam_2 = [0 0 1 1]; qam_6 = [0 1 1 1]; qam_10 = [1 1 1 1]; qam_14 = [1 0 1 1];
        qam_3 = [0 0 1 0]; qam_7 = [0 1 1 0]; qam_11 = [1 1 1 0]; qam_15 = [1 0 1 0];
    
        if (input_1 == qam_0)
            u1 = -3+3i;
        elseif (input_1 == qam_1)
            u1 = -3+1i;
        elseif (input_1 == qam_2)
            u1 = -3-1i;    
        elseif (input_1 == qam_3)
            u1 = -3-3i;
        elseif (input_1 == qam_4)
            u1 = -1+3i;    
        elseif (input_1 == qam_5)
            u1 = -1+1i;    
        elseif (input_1 == qam_6)
            u1 = -1-1i;    
        elseif (input_1 == qam_7)
            u1 = -1-3i;    
        elseif (input_1 == qam_8)
            u1 = 1+3i;    
        elseif (input_1 == qam_9)
            u1 = 1+1i;
        elseif (input_1 == qam_10)
            u1 = 1-1i;    
        elseif (input_1 == qam_11)
            u1 = 1-3i;    
        elseif (input_1 == qam_12)
            u1 = 3+3i;    
        elseif (input_1 == qam_13)
            u1 = 3+1i;    
        elseif (input_1 == qam_14)
            u1 = 3-1i;  
        elseif (input_1 == qam_15)
            u1 = 3-3i; 
        end

        if (input_2 == qam_0)
            u2 = -3+3i;
        elseif (input_2 == qam_1)
            u2 = -3+1i;
        elseif (input_2 == qam_2)
            u2 = -3-1i;    
        elseif (input_2 == qam_3)
            u2 = -3-3i;
        elseif (input_2 == qam_4)
            u2 = -1+3i;    
        elseif (input_2 == qam_5)
            u2 = -1+1i;    
        elseif (input_2 == qam_6)
            u2 = -1-1i;    
        elseif (input_2 == qam_7)
            u2 = -1-3i;    
        elseif (input_2 == qam_8)
            u2 = 1+3i;    
        elseif (input_2 == qam_9)
            u2 = 1+1i;
        elseif (input_2 == qam_10)
            u2 = 1-1i;    
        elseif (input_2 == qam_11)
            u2 = 1-3i;    
        elseif (input_2 == qam_12)
            u2 = 3+3i;    
        elseif (input_2 == qam_13)
            u2 = 3+1i;    
        elseif (input_2 == qam_14)
            u2 = 3-1i; 
        elseif (input_2 == qam_15)
            u2 = 3-3i;
        end
        
%%% 16-QAM signal to transmit  
        signal_qam = [u1; u2];
%%% Channel
        channel = 1/sqrt(2)*(rand(1,2)+1i*rand(1,2));
%%% Alumnati Channel
        alam_channel = [channel(1) channel(2); conj(channel(2)) -conj(channel(1))];
%%% Alumnati Code with channel & noise        
        alam_signal = alam_channel*signal_qam;
        vv = 3.4;  %%% to keep snr=1/No
        alam_noise = vv*(1/sqrt(2))*(randn(2,1) + 1i*randn(2,1));
%%% Signal Received        
        alam_receive_signal = alam_signal + 10^(-Eb_N0_dB(p)/20)*alam_noise;
%%% Alumnati Unitary Channel matrix
        alam_inv_unitary = inv(alam_channel);
%%% Alumnati y~ (y tilda)        
        receive_tilda_y = alam_inv_unitary*alam_receive_signal;
        
%%% ML hard decoding
        received_u1 = receive_tilda_y(1);
        received_u2 = receive_tilda_y(2);
        re_u1 = real(received_u1); im_u1 = imag(received_u1);
        re_u2 = real(received_u2); im_u2 = imag(received_u2);
        op_u1_1 = []; op_u1_2 = []; op_u2_1 = []; op_u2_2 = [];
       
        if (re_u1>=2)
            re_u1 = 3;
            op_u1_1 = [1 0];
        elseif (re_u1>0 && re_u1<2)
            re_u1 = 1;
            op_u1_1 = [1 1];
        elseif (re_u1<=0 && re_u1>-2)
            re_u1 = -1;
            op_u1_1 = [0 1];    
        elseif (re_u1<=-2)
            re_u1 = -3;
            op_u1_1 = [0 0];    
        end
        if (im_u1>=2)
            im_u1 = 3i;
            op_u1_2 = [0 0];
        elseif (im_u1>0 && im_u1<2)
            im_u1 = 1i;
            op_u1_2 = [0 1];
        elseif (im_u1<=0 && im_u1>-2)
            im_u1 = -1i;
            op_u1_2 = [1 1];    
        elseif (im_u1<=-2)
            im_u1 = -3i;
            op_u1_2 = [1 0];    
        end
        op_u1 = [op_u1_1 op_u1_2];

        if (re_u2>=2)
            re_u2 = 3;
            op_u2_1 = [1 0];
        elseif (re_u2>0 && re_u2<2)
            re_u2 = 1;
            op_u2_1 = [1 1];
        elseif (re_u2<=0 && re_u2>-2)
            re_u2 = -1;
            op_u2_1 = [0 1];    
        elseif (re_u2<=-2)
            re_u2 = -3;
            op_u2_1 = [0 0];    
        end
        if (im_u2>=2)
            im_u2 = 3i;
            op_u2_2 = [0 0];
        elseif (im_u2>0 && im_u2<2)
            im_u2 = 1i;
            op_u2_2 = [0 1];
        elseif (im_u2<=0 && im_u2>-2)
            im_u2 = -1i;
            op_u2_2 = [1 1];    
        elseif (im_u2<=-2)
            im_u2 = -3i;
            op_u2_2 = [1 0];    
        end
        op_u2 = [op_u2_1 op_u2_2];
        output = [op_u1 op_u2];

%%% counting the errors 

        % BER
        bererror = bererror + (size(find([input - output]),2)); 
        
        % SER
        if (size(find([op_u1 - input_1]),2) > 0)
            sererror = sererror + 1;
        else
            sererror = sererror;
        end
        if (size(find([op_u2 - input_2]),2) > 0)
            sererror = sererror + 1;
        else
            sererror = sererror;
        end
        
        % CER
        if (size(find([op_u1 - input_1]),2) > 0 || size(find([op_u2 - input_2]),2) > 0)
            cererror = cererror + 1;
        else
            cererror = cererror;
        end            
    end
    cer_err_alum(p) = cererror/itr; % simulated CER
    ser_err_alum(p) = sererror/(2*itr); % simulated SER
    ber_err_alum(p) = bererror/(N*itr); % simulated ber
end
figure;
semilogy(Eb_N0_dB, ber_err_alum,'mo-','LineWidth',2); 
hold on;
semilogy(Eb_N0_dB, ser_err_alum,'bo-','LineWidth',2);
hold on;
semilogy(Eb_N0_dB, cer_err_alum,'go-','LineWidth',2);
axis([Eb_N0_dB(1) Eb_N0_dB(end) 10^-5 0.5]); grid on;
legend('BER Tx=2, Rx=1','SER Tx=2, Rx=1','CER Tx=2, Rx=1');
xlabel('SNR in dB'); ylabel('Probability of error, Pe');
title('Alamouti Code with 16-QAM Constellation by EE18MTECH01005');



