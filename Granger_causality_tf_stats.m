%%  RGS

addpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/fieldtrip')
addpath(genpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/CorticoHippocampal'))
addpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/ADRITOOLS')
addpath('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session')

GC_Bp = GC_Bp_cluster1_rgs;
GC = GC_cluster1_rgs;

I = cellfun(@isempty,GC);
I = not(I);
GC = GC(I);
GC_Bp = GC_Bp(I);


L = cellfun(@length, GC,'UniformOutput', false);
L = cell2mat(L);
GC = GC(L~=1);
GC_Bp = GC_Bp(L~=1);
R = (cellfun(@(equis1) max(abs(hilbert(equis1(2,3001-50:3001+50)))),GC_Bp));
R(R>150) = 0;
[~,r_nl] = sort(R,'descend');
R = R(r_nl);
    
GC = GC(r_nl);
GC = GC(1:1000);

GC_Bp = GC_Bp(r_nl);
GC_Bp = GC_Bp(1:1000);

freqrange = [0:0.5:20];
% freqrange = [20:1:100];
% freqrange = [100:2:300];

p = GC;
q = GC_Bp;

%% Iteration of GC trials RGS
iter = 30;
m = 200;
   grangerspctrm_concat = zeros(2,2,length(freqrange),length(-3:0.01:3),iter);
   for i = 1:iter
       i
       randorder = randperm(length(q));
       temp_q = q(randorder);
       temp_p = p(randorder);
       q_ = temp_q(1:m);
       p_ = temp_p(1:m);
       
       % Compute time frequency gc
fn = 1000;
leng = length(p_);
ro = 3000;
tm = create_timecell(ro,leng);
label = [{'PFC'}; {'HPC'}];

Data.label = label;
Data.time = tm;
Data.trial = p_.';
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data); 
cfg.bsfreq = [249 251];
Data = ft_preprocessing(cfg,Data); 

        [granger_tf] = createauto_timefreq(Data,freqrange);       
        grangerspctrm_concat(:,:,:,:,i) = granger_tf.grangerspctrm;
       
   end
   
 granger_tf.grangerspctrm_RGS_20Hz = grangerspctrm_concat;

%%  Veh

GC_veh_Bp = GC_Bp_cluster1_veh;
GC_veh = GC_cluster1_veh;
I = cellfun(@isempty,GC_veh);
I = not(I);
GC_veh = GC_veh(I);
GC_veh_Bp = GC_veh_Bp(I);

label = [{'PFC'}; {'HPC'}];

L = cellfun(@length, GC_veh,'UniformOutput', false);
L = cell2mat(L);
GC_veh = GC_veh(L~=1);
GC_veh_Bp = GC_veh_Bp(L~=1);
R = (cellfun(@(equis1) max(abs(hilbert(equis1(2,3001-50:3001+50)))),GC_veh_Bp));
R(R>150) = 0;
[~,r_nl] = sort(R,'descend');
R = R(r_nl);
    
GC_veh = GC_veh(r_nl);
GC_veh = GC_veh(1:1000);

GC_veh_Bp = GC_veh_Bp(r_nl);
GC_veh_Bp = GC_veh_Bp(1:1000);

freqrange = [0:0.5:20];
% freqrange = [20:1:100];
% freqrange = [100:2:300];

p = GC_veh;
q = GC_veh_Bp;

%% Iteration of GC trials Veh 
iter = 30;
m = 200;
   grangerspctrm_concat = zeros(2,2,length(freqrange),length(-3:0.01:3),iter);
   for i = 1:iter
       i
       randorder = randperm(length(q));
       temp_q = q(randorder);
       temp_p = p(randorder);
       q_ = temp_q(1:m);
       p_ = temp_p(1:m);
       
% Compute time frequency gc
fn = 1000;
leng = length(p_);
ro = 3000;
tm = create_timecell(ro,leng);
label = [{'PFC'}; {'HPC'}];
Data.label = label;
Data.time = tm;
Data.trial = p_.';
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data); 
cfg.bsfreq = [249 251];
Data = ft_preprocessing(cfg,Data); 

        [granger_tf] = createauto_timefreq(Data,freqrange);
        grangerspctrm_concat(:,:,:,:,i) = granger_tf.grangerspctrm;
       
   end
   
granger_tf.grangerspctrm_Veh_20Hz = grangerspctrm_concat;   
   