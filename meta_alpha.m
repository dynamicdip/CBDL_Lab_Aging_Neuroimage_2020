%Metastability for all 650 subjects in the alpha band 


tic
files = dir('/dataslave2/shyamchand/aamod_meg_maxfilt_00001/CC*');
%files = dir('C:\CamCAN\MiddleElderly\CC*');

meta = zeros(650,1);
fs = 100;
for subs = 1:650;
    try
        dat = ft_read_data(['/dataslave2/shyamchand/aamod_meg_maxfilt_00001/' files(subs).name '/rest/transdef_mf2pt2_rest_raw.fif']);
        %dat = ft_read_data(['C:\CamCAN\MiddleElderly\' files(subs).name '\transdef_mf2pt2_rest_raw.fif']);    
    catch
        continue
    end
    
    dat = dat(3:3:306,:);
   % NChan * Ntime
    dat = resample(dat',100,1000); %  NTime * NChan
    high = 12;
    low = 8;
    filt_sig = ft_preproc_bandpassfilter(dat',100,[high low],[],'fir'); % Channels * Time
    phase_sig = angle(hilbert(filt_sig')); % Ntime * NChannels
    meta_val = metastability(phase_sig');
    meta(subs,1) = meta_val;        
      
    
end
save('meta_alpha_resample','meta')
toc

