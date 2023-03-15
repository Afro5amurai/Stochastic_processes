clc;
clear all ;
close all;
%%
DOF_data = importdata('C:\Users\gilim\Documents\studies\MATLAB\Semester 5\stochastic processes\DOF_measurements.txt');
I_data = readmatrix('C:\Users\gilim\Documents\studies\MATLAB\Semester 5\stochastic processes\I_measurements.txt')';
%%
DOF_data.data=DOF_data.data';
animate_mov_data(DOF_data.data);
% animate_mov_component(DOF_data)
%%
anim_t = linspace(1,120,width(DOF_data.data));
figure(100)
plot(anim_t , DOF_data.data(1:20,:))

%%
%P_neutral
P_neutral = calc_time_avg(DOF_data.data);
dp = DOF_data.data - P_neutral;
dp_dpt = dp*dp';
K_dp = calc_time_avg;
%%
[eig_vec , eig_val] = eig(dp_dpt);
%%

%%
% Define the time axis
rel_rows = 300;
alpha = 0.5;
beta = 0.5;
t_final = 300;
step_size= 1e-2;
standard_deviation = 1;
x_0 = 0;
t_0 = 0;
[x,t] = create_Langevin(x_0,t_0,alpha,t_final,step_size,standard_deviation,0);
len_rel = length(x);
data_matrix =zeros(rel_rows,len_rel) ;
z_arr = zeros(1,len_rel) ;

figure(1)
    hold on
    for i = 1:1:rel_rows
        [x,t,z] = create_Langevin(x_0,t_0,alpha,t_final,step_size,standard_deviation,0);
        data_matrix(i,:) = x;
        z_arr(i,1:end-1) = z;
        plot(t,x)
        title("Realizations of the process")
        xlabel("Time")
        ylabel('Amplitude')
    end
hold off

%%
%Ensemble average calculation
X_ensemble_avg = mean(data_matrix);
figure(2)
plot(t,X_ensemble_avg)
title('Ensemble average')
%%
figure(3)
title('histogram')
histfit(data_matrix(:,end),100)

%%
X_time_avg = calc_time_avg(round(rel_rows/2) , data_matrix);
figure(4)
plot(t, X_time_avg)
%%
X_rms = rms(X_time_avg(:,end));
%%
% 1.1.c
%numeric Autocorrelation calculatuon
% 1
theoretical_R_X_tt = (standard_deviation^2/(2*alpha))*(exp(-alpha*abs(0))-exp(-alpha*(2*t)));
R_X = mean(data_matrix.^2)*step_size;
figure(5)
plot(t,theoretical_R_X_tt)
hold on
plot(t,R_X)
title('Ensemble average for the autocorrelation function of X')
xlim([0 10])
hold off

%%
%2

tau_sample_vec = 0:0.1:10;
ts = 1.5/step_size;
time_mean_X = calc_time_mean(tau_sample_vec , step_size , t , ts , data_matrix , rel_rows);


figure(6)
plot(tau_sample_vec , mean(time_mean_X(:,round(length(time_mean_X)/10):round(length(time_mean_X)/10):length(time_mean_X)),2),"*",'LineWidth',1.5)
hold on
plot(tau_sample_vec , (time_mean_X(:,round(length(time_mean_X)/20):round(length(time_mean_X)/20):round(length(time_mean_X)/2))))
xlim([0 10])
legend("calculation - mean of all realizations" , 'calculation - a few realizations')
hold off
%%
%1.1.D
%numerical sol

[f,peri_ensemble_X , peri_ensemble_Z ,peri_ensemble_X_Z, peri_X ,peri_Z,peri_X_Z] = calc_periodogram(step_size,t,ts,data_matrix,z_arr);

%Ensembel avg   
%%
%analytical sol
figure(7)
semilogx(f , 20*log10(step_size*peri_X(1:4,:)))
title('realizations of PSD X ')
S_X_theory = (standard_deviation/sqrt(step_size))^2*((alpha^2+(2*pi*f).^2)).^-1; 

figure(8)
semilogx(f , 20*log10(step_size*peri_ensemble_X))
hold on
semilogx(f , 20*log10(S_X_theory))
xlim([10^-2 10])
hold off
title('ensemble avg of PSD X')

%%
% 1.1.E
%analytic sol
S_X_Z_theory = (standard_deviation/sqrt(step_size))^2*(alpha+(2j*pi*f)).^-1;
figure(9)
semilogx(f , 20*log10(step_size*abs(peri_X_Z(1:4,:))))
title('realizations of cross PSD X,Z ')
xlim([10^-2 10])

figure(10)
semilogx(f , 20*log10(abs(peri_ensemble_X_Z)))
hold on
semilogx(f , 20*log10(abs(S_X_Z_theory)))
xlim([10^-2 10])
hold off
title('ensemble average of Cross PSD XZ ')
%%
%coherence calc
coherence_X_Z_theo = abs(S_X_Z_theory).^2./(S_X_theory*standard_deviation^2); 
coherence_X_Z = abs(peri_ensemble_X_Z).^2./abs(peri_ensemble_X.*peri_ensemble_Z);
figure(11)
semilogx(f ,coherence_X_Z)
hold on
semilogx(f ,step_size*coherence_X_Z_theo)
hold off
legend('analytical','theoretical')
title('Calculated coherence as a function of frequency')
xlim([10^-2 10])
ylim([0 2])


%%
%------------------------------------------------------------------------
%1.2.A
%1.1.A.1
Y_data_matrix = zeros(height(data_matrix),width(data_matrix));
for j = 1:1:height(data_matrix)
    Y_data_matrix(j,:) = create_Langevin(x_0,t_0,beta,t_final,step_size,standard_deviation,data_matrix(j,:));
end

%%
figure(12)
        plot(t,Y_data_matrix(1:round(rel_rows*0.5),:))
        title("Realizations of the process")
        xlabel("Time")
        ylabel('Amplitude')
figure(13)
title('histogram')
histfit(Y_data_matrix(:,end),100)

%%
%1.1.B.2
Y_ensemble_avg = mean(Y_data_matrix);

figure(14)
plot(t,Y_ensemble_avg)
title('Ensemble average')
Y_time_avg = calc_time_avg(round(rel_rows*0.8) , Y_data_matrix);

figure(15)
plot(t, Y_time_avg)
title('time average as a function of time')
%%
Y_rms = rms(Y_time_avg(:,end));
%%
%1.1.C
% 1
R_Y = mean(Y_data_matrix.^2)*step_size;
figure(16)
plot(t,R_Y)
title('Ensemble average for the autocorrelation function of X')
xlim([0 10])
hold off

% 2
ts = 1.5/step_size;
time_mean_Y = calc_time_mean(tau_sample_vec , step_size , t , ts , Y_data_matrix , rel_rows);

figure(17)
plot(tau_sample_vec , mean(time_mean_Y(:,round(length(time_mean_Y)/10):round(length(time_mean_Y)/10):round(length(time_mean_Y))),2),"-*",'LineWidth',1.5)
hold on
plot(tau_sample_vec , (time_mean_Y(:,round(length(time_mean_Y)/20):round(length(time_mean_Y)/20):round(length(time_mean_Y)/2))))
xlim([0 10])
legend("calculation - mean of all realizations" , 'calculation - a few realizations')
hold off

%1.1.D
[f,peri_ensemble_X,peri_ensemble_Y ,peri_ensemble_X_Y,peri_X,peri_Y,peri_X_Y] = calc_periodogram(step_size,t,ts,data_matrix,Y_data_matrix);
figure(18)
semilogx(f , 20*log10(step_size*peri_Y(1:4,:)))
title('realizations of PSD Y ')
xlim([10^-2 10])

figure(19)
semilogx(f , 20*log10(step_size*peri_ensemble_Y))
xlim([10^-2 10])
title('ensamble average of PSD Y')

%%
%1.1.E
figure(20)
semilogx(f , 20*log10(abs(step_size*peri_X_Y(1:4,:))))
title('realizations of cross PSD X,Y ')
xlim([10^-2 10])

figure(21)
semilogx(f , 20*log10(abs(peri_ensemble_X_Y)))
xlim([10^-2 10])
title('ensemble average of PSD YX')

%%
%coherence
coherence_X_Y = calc_coherence(peri_ensemble_X,peri_ensemble_Y,peri_ensemble_X_Y);
figure(22)
semilogx(f ,coherence_X_Y)
title('Calculated coherence of PSD XY as a function of frequency')
xlim([10^-2 10])
%%
%------------------------------------------------------------------------
%1.2.B
pd=makedist('Lognormal', 'mu', 1, 'sigma', 2);
ln_Z=random(pd,rel_rows,length(t));

%%
%1.1.A
ln_data_matrix = zeros(height(data_matrix),width(data_matrix));
for j = 1:1:height(data_matrix)
    ln_data_matrix(j,:) = create_Langevin(x_0,t_0,beta,t_final,step_size,standard_deviation,ln_Z(j,:));
end
%%
figure(23)
        plot(t,ln_data_matrix(1:round(rel_rows*0.01),:))
        title("Realizations of the process")
        xlabel("Time")
        ylabel('Amplitude')
figure(24)
title('histogram')
histfit(log(ln_data_matrix(:,end)),100)
% fitdist
%%
%1.1.B
ln_ensemble_avg = mean(ln_data_matrix);

figure(25)
plot(t,ln_ensemble_avg)
title('Ensemble average')
ln_time_avg = calc_time_avg(round(rel_rows*0.01) , ln_data_matrix);

figure(26)
plot(t, ln_time_avg)
title('time average as a function of time')
%%
% 1.1.c
%numeric Autocorrelation calculatuon
% 1
tau_sample_vec = 1:0.01:10;
R_X_ln = mean(ln_data_matrix.^2)*step_size;
figure(27)
plot(t,R_X_ln)
title('Ensemble average for the autocorrelation function of X')
% xlim([0 10])

% 2
ts = 1.5/step_size;
time_mean_ln_X = calc_time_mean(tau_sample_vec , step_size , t , ts , ln_data_matrix , rel_rows);


figure(28)
plot(tau_sample_vec , mean(time_mean_ln_X(:,200:200:1000),2),"-*",'LineWidth',1.5)
hold on
plot(tau_sample_vec , (time_mean_ln_X(:,round(length(time_mean_ln_X)/20):round(length(time_mean_ln_X)/20):round(length(time_mean_ln_X)/2))))
% xlim([0 10])
legend("calculation - mean of all realizations" , 'calculation - a few realizations')
hold off
%%
%1.1.D
[f,peri_ensemble_ln_X,peri_ensemble_ln_Z ,peri_ensemble_ln_X_Z,peri_ln_X,peri_ln_Z,peri_ln_X_Z] = calc_periodogram(step_size,t,ts,ln_data_matrix,ln_Z);
figure(29)
semilogx(f , 20*log10(step_size*peri_ln_X(1:4,:)))
title('realizations of PSD X ')
xlim([10^-2 10])

figure(30)
semilogx(f , 20*log10(step_size*peri_ensemble_ln_X))
xlim([10^-2 10])
title('ensamble average of PSD X')
%%
%1.1.E
figure(31)
semilogx(f , 20*log10(abs(step_size*peri_ln_X_Z(1:4,:))))
title('realizations of cross PSD X,Z ')
xlim([10^-2 10])

figure(32)
semilogx(f , 20*log10(abs(peri_ensemble_ln_X_Z)))
xlim([10^-2 10])
title('ensemble average of PSD XZ')

%%
%coherence
coherence_ln_X_Z = calc_coherence(peri_ensemble_ln_X,peri_ensemble_ln_Z,peri_ensemble_ln_X_Z);
figure(33)
semilogx(f ,coherence_ln_X_Z)
title('Calculated coherence of PSD XZ as a function of frequency')
xlim([10^-2 10])
%%
%part 2


%%
%this function solves the langevin equation.
%it takes initial conditions , alpha ,t final, step size , standard
%deviation and an input, for white gaussian noise with mean of zero and std
%of the input standard deviation enter zero.
function[x,t,z] =create_Langevin(x_0,t_0,alpha,t_fin,del_t,sig,input)

    t = t_0:del_t:t_fin;
    sig_d = sig/sqrt(del_t);

    x = zeros(1,length(t));
    z =sig_d*sqrt(del_t)*randn(1,length(t)-1);
    x(1) = x_0;
    t(1) = t_0;    

    if mean(input) ~=0
        for i = 1:length(t)-1
            x(i+1) = x(i) * (1-alpha*del_t) + del_t*input(i);
        end
    else
        for i = 1:length(t)-1
            x(i+1) = x(i) * (1-alpha*del_t) + z(i);
        end
    end
end
function[time_avg] = calc_time_avg(time_avg_samples , data_matrix)
    p = zeros(time_avg_samples, width(data_matrix));
    time_avg = zeros(1,time_avg_samples);
    for i = 1:1:time_avg_samples
        p(i,:) = cumsum(data_matrix(i,:));
        for j = 1:1:length(p)
            time_avg(i,j) = p(i,j)/j;
        end
    end
end
function[f,peri_ensemble_data_1 ,peri_ensemble_data_2 , peri_ensemble_data_1_2 , peri_data_1,peri_data_2,peri_data_1_2] = calc_periodogram(step_size,t,ts,data_1,data_2)
    data_1_after_noise = data_1(:,ts:end);
    transformed_data_1_arr = fft(data_1_after_noise');
    transformed_data_1_arr = transformed_data_1_arr';

    K = 1:1:length(t(ts:end));
    R = 1/step_size;
    N = length(t(ts:end));
    f = (R*K)/N;

    % preallocation 
    peri_data_1 = (1/width(transformed_data_1_arr)) * abs(transformed_data_1_arr).^2;
    peri_ensemble_data_1 = mean(peri_data_1);

    if mean(mean(data_2),2) ~= 0
        data_2_after_noise = data_2(:,ts:end);
        transformed_data_2_arr = fft(data_2_after_noise');
        transformed_data_2_arr = transformed_data_2_arr';
        
        peri_data_1_2 = (1/width(data_1_after_noise)) * (transformed_data_1_arr.*conj(transformed_data_2_arr));
        peri_data_2 = (1/width(transformed_data_2_arr)) * abs(transformed_data_2_arr).^2;
        peri_ensemble_data_1_2 = mean(peri_data_1_2);
        peri_ensemble_data_2 = mean(peri_data_2);
    else
        peri_ensemble_data_1_2 = 0;
        peri_ensemble_data_2 = 0;
        peri_data_1_2 = 0;
        peri_data_2 = 0;
    end
end
function[time_mean] = calc_time_mean(tau_samples,step_size ,t,ts,data_matrix ,realizations)
    tau_vec_normlized = round(tau_samples/step_size);
    time_mean = zeros(length(tau_vec_normlized),length(t));
    for i = 1:1:length(tau_vec_normlized)
        for j = ts:length(t)
            if(j+tau_vec_normlized(i)<width(data_matrix))
                time_mean(i,j) = step_size* mean(data_matrix(1:realizations,j).*data_matrix(1:realizations,j+tau_vec_normlized(i)));
            else
                continue
            end
        end
    end


end
function[coherence_1_2] = calc_coherence(S_1,S_2,S_1_2)
        coherence_1_2 = abs(S_1_2).^2./abs(S_1.*conj(S_2));
end










