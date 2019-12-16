% PRABHAT KUMAR RAI --- EE18MTECH01005
% REPETITIVE PAM & ROTATED PSK

clear all; close all; clc;

snr = 0 : 5 : 30;   % multiple SNR values in dB's

% Generating Input, Noise, Channel & Demodulating the same for different SNR
for p = 1 : length(snr)
    nerror_re_m1 = 0; nerror_re_m2 = 0;
    nerror_ro_m1 = 0; nerror_ro_m2 = 0;
    iteration = 100000; % number of iterations
    
    for q = 1 : iteration
        N = 2; % number of bits or symbols
        input_m1 = rand(1, N)>0.5; % generating 0 or 1 of 2 Dimensional
        input_m2 = rand(1, N)>0.5; % generating 0 or 1 of 2 Dimensional
        input = [input_m1; input_m2];
        
        % Repetitive PAM
        amp = 1/sqrt(10);
        a_1 = [-3*amp; -3*amp] ; b_1 = [-1*amp; -1*amp]; 
        d_1 = [1*amp; 1*amp]; c_1 = [3*amp; 3*amp];
        
        % Rotated PSK
        theta = 0.5*atan(2);
        a_2 = [cos(theta + pi/4); sin(theta + pi/4)];
        b_2 = [cos(theta + 3*pi/4); sin(theta + 3*pi/4)];
        d_2 = [cos(theta - 3*pi/4); sin(theta - 3*pi/4)];
        c_2 = [cos(theta - pi/4); sin(theta - pi/4)];
        
        repetiv_pam_input = []; rotated_psk_input = [];
        
        % Modulation with Rotated PAM and Repetitve PAM
        for pp = 1 : N
            if (input_m1(pp) == 0 && input_m2(pp) == 0)
                repetiv_pam_input = [repetiv_pam_input a_1];
                rotated_psk_input = [rotated_psk_input a_2];
            elseif(input_m1(pp) == 0 && input_m2(pp) == 1) 
                repetiv_pam_input = [repetiv_pam_input b_1];   
                rotated_psk_input = [rotated_psk_input b_2];
            elseif(input_m1(pp) == 1 && input_m2(pp) == 0)      
                repetiv_pam_input = [repetiv_pam_input c_1];      
                rotated_psk_input = [rotated_psk_input c_2];
            elseif(input_m1(pp) == 1 && input_m2(pp) == 1)          
                repetiv_pam_input = [repetiv_pam_input d_1];     
                rotated_psk_input = [rotated_psk_input d_2];
            end
        end
        
        % Channel generation
        h1 = (1/sqrt(2))*[randn(1,1) + 1i*randn(1,1)];
        h2 = (1/sqrt(2))*[randn(1,1) + 1i*randn(1,1)];
        h = [ h1, 0; 0, h2 ];     
        
      % |y1| = | h1 0 |*|repetive_pam_input(1)| + complex| No 0 |
      % |y2|   | 0 h2 | |repetive_pam_input(2)|          | 0 No |
        
        % Channel multiplication
%         repetiv_pam_channel = h*repetive_pam_input;
        repetiv_pam_channel = [h1*repetiv_pam_input(1:2:end); h2*repetiv_pam_input(2:2:end)];
        rotated_psk_channel = [h1*rotated_psk_input(1:2:end); h2*rotated_psk_input(2:2:end)];
        
        % Noise Generation
        noise = (1/sqrt(2))*[randn(2,length(repetiv_pam_channel))+1i*randn(2,length(repetiv_pam_channel))];
        
        % Noise addition
        repetiv_pam_corrup = repetiv_pam_channel + 10^(-snr(p)/20)*noise/sqrt(2);
        rotated_psk_corrup = rotated_psk_channel + 10^(-snr(p)/20)*noise/sqrt(2);
        
        % Equalizing
        repetiv_pam_corrup2 = [repetiv_pam_corrup(1:2:end)*conj(h1)/abs(h1)^2; repetiv_pam_corrup(2:2:end)*conj(h2)/abs(h2)^2];
        rotated_psk_corrup2 = [rotated_psk_corrup(1:2:end)*conj(h1)/abs(h1)^2; rotated_psk_corrup(2:2:end)*conj(h2)/abs(h2)^2];
        
        % Demodulating the signal
        repetiv_pam_corrup2 = real(repetiv_pam_corrup2);
        rotated_psk_corrup2 = real(rotated_psk_corrup2);
        output_re_m1 = []; output_re_m2 = []; output_ro_m1 = []; output_ro_m2 = [];
        
        for pp = 1 : 2 : 2*length(repetiv_pam_corrup2)   
            % Repetitive PAM demodulation via minimum distance decoding rule
            rpd_0 = sqrt((a_1(1)-repetiv_pam_corrup2(pp))^2+(a_1(2)-repetiv_pam_corrup2(pp+1))^2);
            rpd_1 = sqrt((b_1(1)-repetiv_pam_corrup2(pp))^2+(b_1(2)-repetiv_pam_corrup2(pp+1))^2);
            rpd_2 = sqrt((c_1(1)-repetiv_pam_corrup2(pp))^2+(c_1(2)-repetiv_pam_corrup2(pp+1))^2);
            rpd_3 = sqrt((d_1(1)-repetiv_pam_corrup2(pp))^2+(d_1(2)-repetiv_pam_corrup2(pp+1))^2);
            if(rpd_0<rpd_1 && rpd_0<rpd_2 && rpd_0<rpd_3)
                output_re_m1 = [output_re_m1 0];
                output_re_m2 = [output_re_m2 0];
            elseif(rpd_1<rpd_0 && rpd_1<rpd_2 && rpd_1<rpd_3)                
                output_re_m1 = [output_re_m1 0];                
                output_re_m2 = [output_re_m2 1];            
            elseif(rpd_2<rpd_0 && rpd_2<rpd_1 && rpd_2<rpd_3)               
                output_re_m1 = [output_re_m1 1];                
                output_re_m2 = [output_re_m2 0];            
            elseif(rpd_3<rpd_0 && rpd_3<rpd_1 && rpd_3<rpd_2)                
                output_re_m1 = [output_re_m1 1];               
                output_re_m2 = [output_re_m2 1];            
            end                
            % Rotated PSK demdulation via minimum distance decoding rule           
            rrd_0 = sqrt((a_2(1)-rotated_psk_corrup2(pp))^2+(a_2(2)-rotated_psk_corrup2(pp+1))^2);
            rrd_1 = sqrt((b_2(1)-rotated_psk_corrup2(pp))^2+(b_2(2)-rotated_psk_corrup2(pp+1))^2);
            rrd_2 = sqrt((c_2(1)-rotated_psk_corrup2(pp))^2+(c_2(2)-rotated_psk_corrup2(pp+1))^2);
            rrd_3 = sqrt((d_2(1)-rotated_psk_corrup2(pp))^2+(d_2(2)-rotated_psk_corrup2(pp+1))^2);           
            if(rrd_0<rrd_1 && rrd_0<rrd_2 && rrd_0<rrd_3)                
                output_ro_m1 = [output_ro_m1 0];                
                output_ro_m2 = [output_ro_m2 0];            
            elseif(rrd_1<rrd_0 && rrd_1<rrd_2 && rrd_1<rrd_3)
                output_ro_m1 = [output_ro_m1 0];            
                output_ro_m2 = [output_ro_m2 1];            
            elseif(rrd_2<rrd_0 && rrd_2<rrd_1 && rrd_2<rrd_3)        
                output_ro_m1 = [output_ro_m1 1];            
                output_ro_m2 = [output_ro_m2 0];
            elseif(rrd_3<rrd_0 && rrd_3<rrd_1 && rrd_3<rrd_2)              
                output_ro_m1 = [output_ro_m1 1];               
                output_ro_m2 = [output_ro_m2 1];            
            end        
        end
        nerror_re_m1 = nerror_re_m1 + (size(find(input_m1 - output_re_m1),2)); 
        nerror_re_m2 = nerror_re_m2 + (size(find(input_m2 - output_re_m2),2));
        nerror_ro_m1 = nerror_ro_m1 + (size(find(input_m1 - output_ro_m1),2)); 
        nerror_ro_m2 = nerror_ro_m2 + (size(find(input_m2 - output_ro_m2),2)); 
    end
    rep_to_err_m1(p) = nerror_re_m1/(N*iteration);
    rep_to_err_m2(p) = nerror_re_m2/(N*iteration);
    rep_total(p) = (rep_to_err_m1(p)+rep_to_err_m2(p))/2;
    rot_to_err_m1(p) = nerror_ro_m1/(N*iteration);
    rot_to_err_m2(p) = nerror_ro_m2/(N*iteration);
    rot_total(p) = (rot_to_err_m1(p)+rot_to_err_m2(p))/2;
end
% semilogy(EnB, rot_total, 'b*-', 'linewidth', 1.5); hold on;
semilogy(snr, rep_to_err_m1, 'g*-', 'linewidth', 1.5); hold on;
semilogy(snr, rot_to_err_m1, 'r*-', 'linewidth', 1.5); hold on;
axis([0 30 10^-6 0.5]); grid on; legend('Repitive PAM','Rotated PSK'); xlabel('SNR in dB');
ylabel('Probability of error, Pe'); title('Repetitive PAM & Rotated PSK by EE18MTECH01005');
