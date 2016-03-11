
%% ScriptECoGSpectrumModel
%
% Program generates two time-series and powerspectra (PSD) with two different
% settings to simmulate age related changes in the shape of the PSD
%
% For usage, please cite:
%
% Miller KJ, Sorensen LB, Ojemann JG,den Nijs M, (2009). Power-Law Scaling
% in the Brain Surface Electric Potential. PLoS Comput Biol. 5(12):e1000609 
% doi: 10.1371/journal.pcbi.1000609
%
% and
%
% Hermes D. (2015). Levels of synchrony and noise are not the only factors
% that can explain age-related changes in the shape of the power spectrum:
% the effect of changing time-scales of neuronal signaling.
%
%
%     Copyright (C) 2015  D. Hermes & K.J. Miller, Dept. of Psychology New York University and Dept. of Psychology and Dept. of Neurosurgery Stanford University
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%

clear all

%% input parameters

% timesteps, and total time
srate=10000; %sampling rate 
dt=1/srate; %timestep value
t_tot=10; %total time (units of seconds)

nr_sims = 4; % simulation with two sets of parameters

% synaptic stuff, timescale, max firing rate, before modulation
tau_set=[1/(2*pi*75) 1/(2*pi*75) 1/(2*pi*40) 1/(2*pi*40)]; %synaptic decay timescale (units seconds)
alp_set=[1/.05 1/.2 1/.05 1/.2]; % 1/"relaxation time" of broadband across membrane - leakage timescale - units (1/s) - the bullshit part of the model is that this is done prior to the influence of the rhythm
frate=10; %maximum pre-synaptic firing rate (spikes/second) - this isn't appropriately normalized below yet
plot_colors = {'r','m','c','b'};

% number of elements in simulation
num_neurons=6; % number of neurons in simulation
psn_num=6000; %number of presynaptic inputs (roughly 6,000-10,000 for a pyramidal neuron)

%% create pwrlaw input spike pattern, psc timecourse + dendritic potential

%-- loop through multiple simulation settings
data=zeros(t_tot*srate,nr_sims);
tic
figure('Position',[0 0 150 150]), hold on
hold on

% create a spike pattern that is the same across the simulations, then they
% can be directly compared
aa_total = zeros(t_tot*srate,num_neurons);
for nn=1:num_neurons
    toc
    disp(['neuron number ' int2str(nn)])
    %-- convolved spike pattern
    aa=zeros(t_tot*srate,1); % initialize
    for k=1:psn_num % cycle through synapses if using different weights
        a=(rand(t_tot*srate,1))<=(frate/srate); %spike times - think of "a" as the trace at a single synapse
        a=a.*(rand-.5)*2; %arbitrary sign and magnitude on -1 to 1 for subsequent PSC   
        aa=aa+a;
    end
    aa_total(:,nn) = aa;
end

for s = 1:nr_sims
    disp(['starting sim number ' int2str(s)])
    tau_s = tau_set(s);
    alp = alp_set(s);
    
    %-- post-synaptic current shape
    t=0:floor(.05*srate); %50 ms window for PSC
    %create post-synaptic current trace - kjm_pwrlaw_2009 paper way - also see Linden, et. al. 2010
    f=(t.^.13).*exp(1-t/(tau_s*srate)); 
    f=f/sum(f); 
    % plot the post-synaptic corrent trace
    plot(t,f,'Color',plot_colors{s})
    xlim([0 200])
    for nn=1:num_neurons
        toc
        disp(['neuron number ' int2str(nn)])
        %-- convolved spike pattern
        aa=aa_total(:,nn);
        dend_pot=conv(aa,f); %convolve spike arrival  times with PSP shape - superimpose all synapses (i.e. passive dendrite approx.) - same as w/i loop mathematically
        dend_pot((length(a)+1):end)=[];
        
        %-- dendritic temporal integration and leakage, lazy way - Vm is broadband membrane potential
        Vm=0*dend_pot;
        for tt=2:(t_tot*srate)
            Vm(tt)=(1-dt*alp)*Vm(tt-1)+dt*dend_pot(tt-1); %note - this can be off if not careful b/c dend_pot can dominate depending on ratio of std(dend_pot) to alp
        end
        
        data(:,s) = data(:,s) + Vm;
    end %loop through neurons

end % loop through settings

%% calculate power spectra

% [~,f] = psd(data(:,s),srate,srate,srate);
[~,f] = pwelch(data(:,s),srate,srate/2,srate,srate);
spectra_data = zeros(length(f),size(data,2));

% power spectrum
for s = 1:nr_sims
    [spectra_data(:,s),f] = pwelch(data(:,s),srate,srate/2,srate,srate);
%     [spectra_data(:,s),f] = psd(data(:,s),srate,srate,srate);
end


%% make the figure with 4 simulations:

figure('Position',[0 0 200 700])
subplot(3,1,1),hold on
for s = 1:nr_sims
    spect_plot = spectra_data(:,s);
    plot(f,spect_plot,'Color',plot_colors{s})
end
set(gca, 'XScale', 'log', 'YScale', 'log')
set(gca, 'XTick',[1 10 100])
title('loglog')

% set a line at the knees
for s = 1:nr_sims
    f_tau = 1/(2*pi*tau_set(s));
    f_alp = 1/(2*pi*alp_set(s));
    plot([f_tau f_tau],[spect_plot(1) spect_plot(300)],':','Color',plot_colors{s})
    plot([f_alp f_alp],[spect_plot(1) spect_plot(300)],':','Color',plot_colors{s})
end
xlim([1 300])

subplot(3,1,2),hold on
for s = 1:nr_sims
    spect_plot = spectra_data(:,s);
    plot(f,spect_plot,'Color',plot_colors{s})
end
set(gca, 'YScale', 'log')
title('semilog')
xlim([0 300])
set(gca,'XTick',[0 20 50:50:300])
subplot(3,1,3),hold on
for s = 1:nr_sims
    spect_plot = spectra_data(:,s);
    plot(f,spect_plot,'Color',plot_colors{s})
end
set(gca, 'YScale', 'log')
title('semilog')
xlim([0 20])

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/MillerModel_Spectra'])
print('-depsc','-r300',['./figures/MillerModel_Spectra'])

%% make the figure with 2 simulations for the Journal Club

figure('Position',[0 0 200 700])
subplot(3,1,1),hold on
for s = [1 4]%1:nr_sims
    spect_plot = spectra_data(:,s);
    plot(f,spect_plot,'Color',plot_colors{s})

    idx =  f>100 & f<400;
    x = log(f(idx));
    y = log(spect_plot(idx));
    mdl = fitlm(x, y);    
    pwr = mdl.Coefficients.Estimate(2);
    disp(pwr)

end
set(gca, 'XScale', 'log', 'YScale', 'log')
set(gca, 'XTick',[1 10 100])
title('loglog')




% set a line at the knees
for s = [1 4]%1:nr_sims
    f_tau = 1/(2*pi*tau_set(s));
    f_alp = 1/(2*pi*alp_set(s));
    plot([f_tau f_tau],[spect_plot(1) spect_plot(300)],':','Color',plot_colors{s})
    plot([f_alp f_alp],[spect_plot(1) spect_plot(300)],':','Color',plot_colors{s})
end
xlim([1 300])

subplot(3,1,2),hold on
for s = [1 4]%1:nr_sims
    spect_plot = spectra_data(:,s);
    plot(f,spect_plot,'Color',plot_colors{s})
end
set(gca, 'YScale', 'log')
title('semilog')
xlim([0 300])
set(gca,'XTick',[0 20 50:50:300])
subplot(3,1,3),hold on
for s = [1 4]%1:nr_sims
    spect_plot = spectra_data(:,s);
    plot(f,spect_plot,'Color',plot_colors{s})
end
set(gca, 'YScale', 'log')
title('semilog')
xlim([0 20])

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/MillerModel_Spectra_JC'])
print('-depsc','-r300',['./figures/MillerModel_Spectra_JC'])

