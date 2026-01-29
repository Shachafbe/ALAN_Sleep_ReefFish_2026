% routine to run perliminary analysis of the fish movies 

% define source folder where data files (.csv) and movie files (.avi) are
% at. IMPORTANT: movie files and data files should have the same exact
% names, with only different extentions.

% define colormaps for lines 

Cmap = LineMap;

% folder where .csv and .mp4 files are places
% change to your folder
Folder = ''; 

% put here the folder where you want files to be saved to (right now only
% .fig files are saved)
save_folder = '';

% important decide if to save files and make a movie as default (I usually
% leave it as 0, and change it only if I need to).
SAVE = 0;
MAKEMOVIE = 0; % decide if to make a movie or not



cd(Folder); % change to folder location


% put in file names (without extensions):
%filenames_light = {'LN1DeepCut_resnet50_lightdarkNov13shuffle1_650000',...
   %             'LN2DeepCut_~1930-0500',...
  %              'LN3DeepCut_resnet50_lightdarkNov13shuffle1_650000'};
%filenames_dark = {'DN1DeepCut_~1930-0540',...
 %   'DN2DeepCut_resnet50_lightdarkNov13shuffle1_350000'};

% enter recording time in sec. 
% THIS IS IMPORTANT - it will interpulate the 12 fps movies to 13 fps
% according to this value
 recording_time = 180; % in sec

low_fr = 12; % low and high frame rates
tt_low = 1/low_fr:1/low_fr:180;
high_fr = 13;
tt_high = 1/high_fr:1/high_fr:180;

% put in file names (without extensions):
filenames_light = dir('WhiteShortT*.csv'); % all files starting with LN and a .csv extension
filenames_light = {filenames_light.name};
filename_long = dir('WhiteLongT*.csv');
filename_long = {filename_long.name};
filenames_light = [filenames_light, filename_long];
filenames_dark = dir('WhiteBeforeT*.csv'); % all files starting with DN and a .csv extension
filenames_dark = {filenames_dark.name};



orderl = [1,2]; % days of the experiment
orderd = [1,2]; 

numl = length(filenames_light); % number of light experiments
numd = length(filenames_dark);  % number of dark experiments

% lood data files and create movie objects:
coorl = cell(1,numl); 
coord = cell(1,numd);

% create also links to video files
vid_l = cell(1,numl); 
vid_d = cell(1,numd);

Nf = 8; % number of fish (no including the coordinated used for measuring in the files)
Num_dots = 6; % number of markers on a fish
% Tl = size(coorl,1);
% Td = size(coord,1);
output_per_dot = 3; % number of outputs from DLC per dot (x.y.q) just in case this ever changes

% loop to load data files and video objects (takes some time)
for i = 1:numl
    coorl{i} = csvread([Folder,filenames_light{i}],3,1);
    coorl{i} = rearrangeCoorData(coorl{i},'num_fish',Nf,'data_per_point',output_per_dot,...
        'points_per_fish',Num_dots);
   
    if exist([Folder,filenames_light{i}(1:end-4),'_id_labeled.mp4'],'file')
        vid_l{i} = VideoReader([Folder,filenames_light{i}(1:end-4),'_id_labeled.mp4']);
    end
    disp(['Finished light file #',num2str(i)]);
end


for i = 1:numd
    coord{i} = csvread([Folder,filenames_dark{i}],3,1);
    coord{i} = rearrangeCoorData(coord{i},'num_fish',Nf,'data_per_point',output_per_dot,...
        'points_per_fish',Num_dots);
    if exist([Folder,filenames_dark{i}(1:end-4),'_id_labeled.mp4'],'file')
        vid_d{i} = VideoReader([Folder,filenames_dark{i}(1:end-4),'_id_labeled.mp4']);
    end
    disp(['Finished dark file #',num2str(i)]);
end
    
% titles for plotting
Titles_light = {'LongT1C1','LongT5C1'};
Titles_dark = {'BeforeT1C2','BeforeT5C2'};

ref_dot_to_use = 4; % which dot to choose
%% Try to clean up the data a little bit:

% extract head coordinates and quality ratings from the data matrices:

% light data
xhl = cell(1,numl); % head x coor of ligth files
yhl = cell(1,numl); % y coor
qhl = cell(1,numl); % quality of tracking

Shl = cell(1,numl); % variable for calculated speed

for i = 1:numl
    xhl{i} = coorl{i}(:,:,ref_dot_to_use,1); % extract only relevant columns
    yhl{i} = coorl{i}(:,:,ref_dot_to_use,2);
    qhl{i} = coorl{i}(:,:,ref_dot_to_use,3);
    
    % interpulate if needed
    if length(xhl{i}) < length(tt_high)
        xhl{i} = interp1(tt_low,xhl{i},tt_high); % interpulate for high frame rate;
        yhl{i} = interp1(tt_low,yhl{i},tt_high); % interpulate for high frame rate;
        qhl{i} = interp1(tt_low,qhl{i},tt_high); % interpulate for high frame rate;
    end
    
    vx = balancedDerivative(xhl{i});
    vy = balancedDerivative(yhl{i});
    
    Shl{i} = calculateNorm(cat(3,vx,vy));
end


% dark data
xhd = cell(1,numd); % head x coor of ligth files
yhd = cell(1,numd); % y coor
qhd = cell(1,numd); % quality of tracking

Shd = cell(1,numl); % variable for calculated speed

for i = 1:numd
    xhd{i} = coord{i}(:,:,ref_dot_to_use,1); % extract only relevant columns
    yhd{i} = coord{i}(:,:,ref_dot_to_use,2);
    qhd{i} = coord{i}(:,:,ref_dot_to_use,3);
    
     % interpulate if needed
    if length(xhd{i}) < length(tt_high)
        xhd{i} = interp1(tt_low,xhd{i},tt_high); % interpulate for high frame rate;
        yhd{i} = interp1(tt_low,yhd{i},tt_high); % interpulate for high frame rate;
        qhd{i} = interp1(tt_low,qhd{i},tt_high); % interpulate for high frame rate;
    end

    vx = balancedDerivative(xhd{i});
    vy = balancedDerivative(yhd{i});
    
    Shd{i} = calculateNorm(cat(3,vx,vy));
end


% decide speed and waulity thresholds to interpulate over:
%(thesse can be manually changed or later decided upon in a more
%informative way/smarter way)
speed_th = 50; % above this speed a point is considered spurious
quality_th = 0.8; % below this quality a point is considered spurious
filt_win = 5; % number of frames to smooth over 

disp([' removing points with speed > ',num2str(speed_th),...
    ' and quality < ',num2str(quality_th)]);

% clean up the data by interpulating over spurious points:

% variables foe interpulated data
xcl = cell(1,numl); % x for cleaned up light data
ycl = cell(1,numl); % y
Scl = cell(1,numl); % speed
time_corrected_l = cell(Nf,numl); % time points corrected

xcd = cell(1,numd); % x for clearned up dark data
ycd = cell(1,numd); % y
Scd = cell(1,numd); % speed
time_corrected_d = cell(Nf,numd); % time points corrected

% loop over experiments and then over fish:

% for light experiments:
figure;
p = 1; ax = [];
num_cols = max([numl,numd]);
for j = 1:numl % all experiments
    for i = 1:Nf
        
        x = xhl{j}(:,i); y = yhl{j}(:,i);
        
        % find suspicious points to correct
        jj =  (Shl{j}(:,i) > speed_th | qhl{j}(:,i) < quality_th...
            | xhl{j}(:,i)==0 | yhl{j}(:,i) ==0);
        
        time_corrected_l{i,j} = jj;

        jj = find(jj);
        % start interpulating before and after the points:
        temp_qs = false(size(x));
        for pp = -2:2
            ii = jj+pp;
            ii(ii<1 | ii > length(x)) = [];
            temp_qs(ii) = true;
        end
                
        % interpulate over suspicious points:
        x(temp_qs) = interp1(find(~temp_qs),x(~temp_qs),find(temp_qs));
        y(temp_qs) = interp1(find(~temp_qs),y(~temp_qs),find(temp_qs));
        
        % smooth the data 
        xcl{j}(:,i) = smoothdata(x,'sgolay',filt_win);
        ycl{j}(:,i) = smoothdata(y,'sgolay',filt_win);
        
    end
    % calculate speed estimation of the corrected data:
    vx = balancedDerivative(xcl{j});
    vy = balancedDerivative(ycl{j});
    
    Scl{j} = calculateNorm(cat(3,vx,vy));
    
    disp(['Finished correcting light experiment #',num2str(j)]);
    
    % plot amount of time corrected;
    ax(p) = subplot(2,num_cols,j); p = p+1;
    pct = cellfun(@mean,time_corrected_l(:,j));
    bar(pct)
    title(Titles_light{j});
    xlabel('fish'); ylabel('% data corrected');
    box off;
end

% for dark experiments:
for j = 1:numd % all experiments
    for i = 1:Nf
        
        x = xhd{j}(:,i); y = yhd{j}(:,i);
        
        % find suspicious points to correct
        jj =  (Shd{j}(:,i) > speed_th | qhd{j}(:,i) < quality_th...
            | xhd{j}(:,i)==0 | yhd{j}(:,i) ==0);
        
        time_corrected_d{i,j} = jj;

        jj = find(jj);
        % start interpulating before and after the points:
        temp_qs = false(size(x));
        for pp = -2:2
            ii = jj+pp;
            ii(ii<1 | ii > length(x)) = [];
            temp_qs(ii) = true;
        end
        

        
        % interpulate over suspicious points:
        x(temp_qs) = interp1(find(~temp_qs),x(~temp_qs),find(temp_qs));
        y(temp_qs) = interp1(find(~temp_qs),y(~temp_qs),find(temp_qs));
        
        % smooth the data 
        xcd{j}(:,i) = smoothdata(x,'sgolay',filt_win);
        ycd{j}(:,i) = smoothdata(y,'sgolay',filt_win);
        
    end
    
    % calculate speed estimation of the corrected data:
    vx = balancedDerivative(xcd{j});
    vy = balancedDerivative(ycd{j});
    
    Scd{j} = calculateNorm(cat(3,vx,vy));
    disp(['Finished correcting dark experiment #',num2str(j)]);
    % plot amount of time corrected;
    ax(p) = subplot(2,num_cols,j+num_cols); p = p+1;
    pct = cellfun(@mean,time_corrected_d(:,j));
    bar(pct);
    title(Titles_dark{j});
    xlabel('fish'); ylabel('% data corrected');
    box off;
end

linkaxes(ax,'xy');
% plot amount of data we had to correct:
%% plot trajectories and speed distributions:

% choose which groups to plot:
groups_l = [1,2]; % out of 3 groups (i.e. [1,2,3] will plot all 3 light groups, and [] will not plot any)
groups_d = [1,2]; % out of two groups 

% shades for histogram
shades_l = linspace(1,0.2,numl);
shades_d = linspace(1,0.2,numd);

num_rows = max([length(groups_d),length(groups_l)]); % number of columns in subplots
num_cols = 2; % always 2 since we compare dark and light

% th_quality_for_traj_plot = 0.90;
% position map:
% start ploting light night groups (from the left)
f1 = figure();
p = 1;
for l = groups_l % all groups
    subplot(num_rows,num_cols,(p-1)*2+1);
    if ~isempty(vid_l{l})
    img = readFrame(vid_l{l});
    imshow(img); 
    end
    hold on;
    for f = 1:Nf
        % remove any spurious data that might clutter  the plots (just for visualizations):
%         
%         qsl = Scl{l}(:,f) < 20 & qhl{l}(:,f) > th_quality_for_traj_plot;
        jj = time_corrected_l{f,l};
        plot(xcl{l}(~jj,f),ycl{l}(~jj,f),'.','MarkerSize',1);

    end
    title(Titles_light{l});
    p = p+1;
end

% add dark night groups (from the right)
p = 1;
for d = groups_d % all groups
    subplot(num_rows,num_cols,(p-1)*2+2);
    if ~isempty(vid_d{d})
    img = readFrame(vid_d{d});
    imshow(img); 
    end
    hold on;
    for f = 1:Nf
%         qsd = Scd{d}(:,f) < 12 & qhd{d}(:,f) > th_quality_for_traj_plot;
        jj = time_corrected_d{f,d};
        plot(xcd{d}(~jj,f),ycd{d}(~jj,f),'.','MarkerSize',1);

    end
    title(Titles_dark{d});
    p = p+1;
end


if SAVE
    name = ['trajectories_',[Titles_light{groups_l}],...
        [Titles_dark{groups_d}]];
    saveas(gcf,[save_folder,name]);
    disp('saved group trajectories');
end    

% now plot the distributions of speeds:

f2 = figure();

% first light distributions:
p = 1;ax = [];
for l = groups_l
    for f = 1:Nf
        
        ax(f) = subplot(1,Nf,f);
        jj = time_corrected_l{f,l};
        histogram(log(Scl{l}(~jj,f)),-7:0.1:7,'Normalization','Probability',...
            'EdgeColor','none','FaceColor',Cmap(1,:)*shades_l(p)); hold on;
        box off;
        xlabel('log(Speed)');
        title(['Fish ',num2str(f)]);
    end
    p = p+1;
end

% now dark nigth distributions
p = 1;
for d = groups_d
    for f = 1:Nf
        
        subplot(1,Nf,f);
        jj = time_corrected_d{f,d};

        histogram(log(Scd{d}(~jj,f)),-7:0.1:7,'Normalization','Probability',...
            'EdgeColor','none','FaceColor',Cmap(2,:)*shades_d(p)); hold on;
        box off;
        xlabel('log(Speed)');
    end
    p = p+1;
end

% add legend to last subplot
legend([Titles_light(groups_l), Titles_dark(groups_d)]);
legend boxoff;
linkaxes(ax,'xy');
%     subplot(2,3,i);
        %     plot(Vorig); hold on;
        %     plot(V);
        %     plot(V); hold on;
        %     histogram(V);

if SAVE
    name = ['speed_distributions_',[Titles_light{groups_l}],...
        [Titles_dark{groups_d}]];
    saveas(gcf,[save_folder,name]);
    disp('saved speed distributions');
end
%% look at times of 'up' and 'down' states:

% define thresholds for "up" and "down" states:
% th = exp([0.01 0.01 0.01]);
log_th = 0.4; % previously was 0.1
th = exp(ones(1,Nf)*log_th);

% for simplicity we use the same th for all fish in all experiments. but we
% can make a smarter desicion according to the speed distributions of the
% experiments we are comparing

% variables for up down states
up_down_l = cell(1,numl);
up_down_d = cell(1,numd);

% variables for start and end of up down of events (time frames), and for
% there length (end-start+1):
lims_l = cell(1,numl);
length_l = cell(1,numl);

lims_d = cell(1,numd);
length_d = cell(1,numd);

% variables for length of low quality events (control analysis): 
length_q_l = cell(3,1); % for length of non quality events
length_q_d = cell(3,1); % for length of non quality events

% define short events to ignore:
short_events = 1; % in frames

% light experiments:
for l = 1:numl
    % temporary varibles per experiment
    temp_up_down = false(size(Scl{l}));
    temp_length = cell(1,Nf); % for 3 fish
    temp_lims = cell(1,Nf);
    
    for f = 1:Nf % loop over fish
        temp = Scl{l}(:,f); % take speed
        jj = time_corrected_l{f,l};
%         temp(jj) = -1000; % remove points that are interpulated

        up_down = temp > th(f); % find up and down
        [st,nd,temp_up_down(:,f)] = findNonZeroSeq(up_down, short_events); % find non zero sequences (remove short events)
        temp_length{f} = nd-st+1; % calcualte their length
        temp_lims{f} = [st;nd];
        
    end
    % save all fish data in experiment variables:
    up_down_l{l} = temp_up_down;
    length_l{l} = temp_length;
    lims_l{l} = temp_lims;
    
end


% dark experiments:
for d = 1:numd
    % temporary varibles per experiment
    temp_up_down = false(size(Scd{d}));
    temp_length = cell(1,Nf); % for 3 fish
    temp_lims = cell(1,Nf);
    
    for f = 1:Nf % loop over fish
        temp = Scd{d}(:,f); % take speed
        j = time_corrected_d{f,d};
%         temp(jj) = -1000; 
        
        up_down = temp > th(f); % find up and down
        [st,nd,temp_up_down(:,f)] = findNonZeroSeq(up_down, short_events); % find non zero sequences (remove short events)
        temp_length{f} = nd-st+1; % calcualte their length
        temp_lims{f} = [st;nd];
        
    end
    % save all fish data in experiment variables:
    up_down_d{d} = temp_up_down;
    length_d{d} = temp_length;
    lims_d{d} = temp_lims;
end
%% plot the distributions of up and down states

figure();

max_ylim = 50; % choose the max values of y axis to plot (to make the grpah look good)
max_xlim = 100; % maximum numbe of bins to show in x axis
% redefine groups to plot:
% choose which groups to plot:
groups_l = [1,2]; % out of 3 groups (i.e. [1,2,3] will plot all 3 light groups, and [] will not plot any)
groups_d = [1,2]; % out of two groups 

% shades for histogram
shades_l = linspace(1,0.2,numl);
shades_d = linspace(1,0.2,numd);
ax = [];
p = 1;
for l = groups_l
    for f = 1:Nf
        
        ax(f) = subplot(1,Nf,f);
        h = histogram(length_l{l}{f},0:1:max_xlim,...
            'EdgeColor','none','FaceColor',Cmap(1,:)*shades_l(p));
        hold on;
        box off;
        if f==round(Nf/2)
        xlabel('Length of activity events [frames]');
        end
    end
end


p = 1;
for d = groups_d
    for f = 1:Nf
        
        subplot(1,Nf,f);
        histogram(length_d{d}{f},0:1:max_xlim,...
            'EdgeColor','none','FaceColor',Cmap(2,:)*shades_d(p));
        hold on;
        box off;
    end
end

% add legend to last subplot
legend([Titles_light(groups_l), Titles_dark(groups_d)]);
legend boxoff;

linkaxes(ax,'xy');
%ylim([0 max_ylim]);

if SAVE
    name = ['active_in_active_dist_',[Titles_light{groups_l}],...
        [Titles_dark{groups_d}]];
    saveas(gcf,[save_folder,name]);
    disp('saved active and inactive distributions');
end


%% explore different sleep th values on first dark and light nights:

sleep_ths = [1:1:10]; % choose various sleep th's in seconds
speed_th = mean(th);


T_min = (1:length(Scl{1}))./Fs/60;
T_min = floor(T_min(end));



ii_l = [1,2]; % index of light movie
ii_d = [1,2]; % index of dark movie

all_total_l = zeros(Nf,length(sleep_ths),length(ii_l));
all_total_d = zeros(Nf,length(sleep_ths),length(ii_d));

all_sleep_pct_l = zeros(T_min,Nf,length(sleep_ths),length(ii_l)); 
all_sleep_pct_d = zeros(T_min,Nf,length(sleep_ths),length(ii_d)); 
for i = 1:length(sleep_ths)



    curr_th = sleep_ths(i);

    for j = 1:length(ii_l)
        [all_sleep_pct_l(:,:,i,j),all_total_l(:,i,j)] = measureSleep(Scl{ii_l(j)},time_corrected_l(:,ii_l(j)),...
            'speed_th',speed_th,'sleep_th',curr_th,'fps',Fs,...
            'sleep_per_time','min');
    end

    for j = 1:length(ii_d)
        [all_sleep_pct_d(:,:,i,j),all_total_d(:,i,j)] = measureSleep(Scd{ii_d(j)},time_corrected_d(:,ii_d(j)),...
            'speed_th',speed_th,'sleep_th',curr_th,'fps',Fs,...
            'sleep_per_time','min');
    end





    i
end

% now plot pct of time sleeping for ever th and for both days:
figure;
fracs = [1.5 1 0.75 0.5];
for j = 1:size(all_total_l,3)
    curr_color = Cmap(1,:)*fracs(j);
    curr_color(curr_color>1) = 1;
h1(j) = boundedLine(sleep_ths,mean(all_total_l(:,:,j)),stdE(all_total_l(:,:,j)),...
    'cmap',curr_color,'alpha');
hold on;
end
for j = 1:size(all_total_d,3)
    curr_color = Cmap(2,:)*fracs(j);
    curr_color(curr_color>1) = 1;
h2(j) = boundedLine(sleep_ths,mean(all_total_d(:,:,j)),stdE(all_total_d(:,:,j)),...
    'cmap',curr_color,'alpha');
hold on;
end
xlabel('Sleep threshold [s]');
ylabel('% of time sleeping');
legend([h1,h2],{Titles_light{ii_l},Titles_dark{ii_d}}); legend box off;

%% also just plot a box plot for the total sleep (as it doesn't male since to look 'per minute'

sleep_th = 7; % choose a single th to plot full night data

ii = find(sleep_ths<=sleep_th,1,'last');
actual_th = sleep_ths(ii);

temp_l = all_sleep_pct_l(:,:,ii);
temp_d = all_sleep_pct_d(:,:,ii);


temp_l = squeeze(all_total_l(:,ii,:));
temp_d = squeeze(all_total_d(:,ii,:));

figure;
h = boxplot([temp_l,temp_d]);
hold on;
for j = 1:size(temp_l,2)
plot(ones(1,length(temp_l))*j,temp_l(:,j),'.','color',Cmap(j,:))
end
for j = 1:size(temp_d,2)
plot(j+size(temp_l,2),(temp_d(:,j)),'.','color',Cmap(j+size(temp_l,2),:));
end




ylabel('% of time sleeping');
set(gca,'XTickLabel',[Titles_light,Titles_dark])

box off;
all_data = [temp_l,temp_d];
all_names = [Titles_light,Titles_dark];
data_tbl = table('Size',size(all_data),'VariableTypes',repmat({'double'},1,size(all_data,2)),...
    'VariableNames',all_names);

for i = 1:size(all_data,2)
    data_tbl(:,i) = table(all_data(:,i));
end

filename = ['Percent of time sleeping sleeph TH is ',num2str(sleep_th) ' seconds.csv'];
writetable(data_tbl,filename)
save(filename(1:end-4),'data_tbl')

%% calculate sleep bout length, wake bout length, # of transitions

Fs = 13; % frames per seconds in movie
sleep_th = 7*Fs; % in frames
connect_th = 1; % below this number of frames - connect sleep sequences 

% time sleeping in seconds/minute 
sleep_rate_l = cell(1,numl);
sleep_rate_d = cell(1,numd);

% mean length of sleeping bouts in seconds
sleep_bout_length_l = cell(1,numl);
sleep_bout_length_d = cell(1,numd);

% mean time of active bouts in seconds
awake_bout_length_l = cell(1,numl);
awake_bout_length_d = cell(1,numd);

% rate of transitions between states in transitions/minute
sleep_rate_transitions_l = cell(1,numl);
sleep_rate_transitions_d = cell(1,numd);


% in % of night
for l = 1:numl
  

    % find # of tranaitions from sleep to wake

    temp = up_down_l{l};
    new_temp = zeros(size(temp))';
    rate_of_transitions = zeros(1,Nf);
    sleep_bout_length = zeros(1,Nf);
    awake_bout_length = zeros(1,Nf);
    sleep_rate = zeros(1,Nf);
    for f = 1
    [st,nd,new_temp(f,:)] = findNonZeroSeq(~temp(:,f),sleep_th,connect_th);
    L = nd-st+1;
    sleep_bout_length(f) = mean(L/Fs);
    sleep_rate(f) = sum(L)/length(temp)*60;

    % find awake bout length (anything that is not sleep in this caese)
    [st,nd] = findNonZeroSeq(~new_temp(f,:));
    L = nd-st+1;
    awake_bout_length(f) = mean(L/Fs);
    % find number of transitions
    diff_new_temp = diff(new_temp(f,:),1,2);
    rate_of_transitions(f) = sum(abs(diff_new_temp)>0)/(length(new_temp)/Fs)*60; % transitions per min
 
    end
    sleep_rate_transitions_l{l} = rate_of_transitions;
    
    sleep_rate_l{l} =sleep_rate;
    
    sleep_bout_length_l{l} = sleep_bout_length;
    awake_bout_length_l{l} = awake_bout_length;
end

% in % of night
for d = 1:numd
  

    % find # of tranaitions from sleep to wake

    temp = up_down_d{d};
    new_temp = zeros(size(temp))';
    rate_of_transitions = zeros(1,Nf);
    sleep_bout_length = zeros(1,Nf);
    awake_bout_length = zeros(1,Nf);
    sleep_rate = zeros(1,Nf);
    for f = 1:Nf
    [st,nd,new_temp(f,:)] = findNonZeroSeq(~temp(:,f),sleep_th,connect_th);
    L = nd-st+1;
    sleep_bout_length(f) = mean(L/Fs);
    sleep_rate(f) = sum(L)/length(temp)*60;

    % find awake bout length (anything that is not sleep in this caese)
    [st,nd] = findNonZeroSeq(~new_temp(f,:));
    L = nd-st+1;
    awake_bout_length(f) = mean(L/Fs);
    % find number of transitions
    diff_new_temp = diff(new_temp(f,:),1,2);
    rate_of_transitions(f) = sum(abs(diff_new_temp)>0)/(length(new_temp)/Fs)*60; % transitions per min
 
    end
    sleep_rate_transitions_d{d} = rate_of_transitions;
    
    sleep_rate_d{d} =sleep_rate;
    
    sleep_bout_length_d{d} = sleep_bout_length;
    awake_bout_length_d{d} = awake_bout_length;
end
%% plot fish heatmaps

% for light groups:
simple_maps = LineMap;
nbins = 50; % number of bins in 2d histogram
smooth_val = 0.5; % smoothing factor for 2D histogram

for l = groups_l
    if ~isempty(vid_l{l})
        img = readFrame(vid_l{l});
    else
        img = [];
    end
    figure('Units','Normalized','Position',[0.0 0.0 1 1]);

    for f = 1:Nf
        tempx = xcl{l}(:,f); % get current coordinates:
        tempy = ycl{l}(:,f);
        
        % remove nans
        tempx = tempx(~isnan(tempx));
        tempy = tempy(~isnan(tempy));

        
        subplot(Nf,2,f*2-1);
        
        if ~isempty(img)
            imshow(img);
        end
        hold on; % draw current image
%         scatter(tempx(1:5:end),tempy(1:5:end),0.1,'MarkerEdgeColor','none',...
%             'MarkerFaceColor',Cmap(f,:),'MarkerFaceAlpha',0.1);
        
        % calculate 2D smooth histogram
        [ctr1,ctr2,Counts] = smoothhist2D([tempx tempy],smooth_val,[nbins nbins],0,'none'); 
%        
        [Nrow,Ncol,~] = size(Counts); % number of rows and columns of histogram       
        % create a color mask same size as histogram
        Color = cat(3,repmat(simple_maps(f,1),Nrow,Ncol),repmat(simple_maps(f,2),Nrow,Ncol),...
            repmat(simple_maps(f,3),Nrow,Ncol));
        % draw the color mask
        hh = imagesc(ctr1,ctr2,Color);colormap(CmapRed);
        % change the transperancy of color mask according to 2d hist such
        % that we will see colors only when heat map is high:
        set(hh,'AlphaData',Counts);
        
        title([Titles_light{l},' fish: ',num2str(f)]);
        % draw only a heat map next to the image for better view
        subplot(Nf,2,f*2); 
        hh = imagesc(ctr1,ctr2,Counts);colormap(CmapRed); hold on;
        set(gca,'Xtick',[],'XtickLabel',[],'Ytick',[],'YtickLabel',[]);
        axis image;
        if ~isempty(img)
        axis([0 size(img,2) 0 size(img,1)]); 
        end
        
    end
    if SAVE
       name = ['Heatmaps_',Titles_light{l}];
       saveas(gcf,[save_folder,name]);
       disp(['saved heatmap for ',Titles_light{l}]);
    end
end

% for dark groups:


for d = groups_d
    if ~isempty(vid_d{d})
        img = readFrame(vid_d{d});
    else
        img = [];
    end
    figure('Units','Normalized','Position',[0.0 0.0 1 1]);
    for f = 1:Nf
        tempx = xcd{d}(:,f); % get current coordinates:
        tempy = ycd{d}(:,f);
        
        tempx = tempx(~isnan(tempx)); % get rid of nans
        tempy = tempy(~isnan(tempy));
        
        subplot(Nf,2,f*2-1);
        
        if ~isempty(img)
            imshow(img); 
        end
        hold on; % draw current image
%         scatter(tempx(1:5:end),tempy(1:5:end),0.1,'MarkerEdgeColor','none',...
%             'MarkerFaceColor',Cmap(f,:),'MarkerFaceAlpha',0.1);
        
        % calculate 2D smooth histogram
        [ctr1,ctr2,Counts] = smoothhist2D([tempx tempy],smooth_val,[nbins nbins],0,'none'); 
%        
        [Nrow,Ncol,~] = size(Counts); % number of rows and columns of histogram       
        % create a color mask same size as histogram
        Color = cat(3,repmat(simple_maps(f,1),Nrow,Ncol),repmat(simple_maps(f,2),Nrow,Ncol),...
            repmat(simple_maps(f,3),Nrow,Ncol));
        % draw the color mask
        hh = imagesc(ctr1,ctr2,Color);colormap(CmapRed);
        % change the transperancy of color mask according to 2d hist such
        % that we will see colors only when heat map is high:
        set(hh,'AlphaData',Counts);
        
        title([Titles_dark{d},' fish: ',num2str(f)]);
        % draw only a heat map next to the image for better view
        subplot(Nf,2,f*2); 
        hh = imagesc(ctr1,ctr2,Counts);colormap(CmapRed); hold on;
        set(gca,'Xtick',[],'XtickLabel',[],'Ytick',[],'YtickLabel',[]);
        axis image;
        if ~isempty(img)
            axis([0 size(img,2) 0 size(img,1)]); 
        end
        
    end
    if SAVE
       name = ['Heatmaps_',Titles_dark{d}];
       saveas(gcf,[save_folder,name]);
       disp(['saved heatmap for ',Titles_dark{d}]);
    end
end
%% plot examples for active and inactive events:

le = 1; % light experiment
de = 1; % dark experiment
f_num = 2;  % fish number

ax = [];
f = figure('Units','Normalized','Position',[0.1 0.3 0.8 0.6]);
ax(1) = subplot(2,1,1);
plot(Scl{le}(:,f_num)); hold on;
backgroundBoxes(gcf,[lims_l{le}{f_num}(1,1:100); lims_l{le}{f_num}(2,1:100)],Cmap(1,:),0.3);
box off; title('Light');
xlabel('Time [frames]'); ylabel('Speed [pix/frames]');

ax(2) = subplot(2,1,2);
plot(Scd{de}(:,f_num),'Color',Cmap(2,:)); hold on
backgroundBoxes(gcf,[lims_d{de}{f_num}(1,1:100); lims_d{de}{f_num}(2,1:100)],Cmap(2,:),0.3);
box off; title('Dark');
xlabel('Time [frames]'); ylabel('Speed [pix/frames]');

linkaxes(ax,'xy');
xlim([0 2000]); ylim([0 50]);

if SAVE
    name = ['example_active_inactive'];
    saveas(gcf,[save_folder,name]);
    disp(['saved example for active and inactive states']);
end
%% make example movie to compare fish speeds with tracking

le = 1; % light experiment
de = 1; % dark experiment
fs = 2;  % fish number
fish_i_for_speed = 2; % fish number for speed plot
% create movie objects
movl = vid_l{le}; 
movd = vid_d{de};


if MAKEMOVIE
        movie_name = ['compare_light_dark_and_speed_fish1_smoothed'];
        obj = VideoWriter([movie_name],'MPEG-4');
        obj.FrameRate = 25;
        obj.Quality = 100;
        open(obj);
end

figure('Units','Normalized','Position',[0.0 0.0 1 1]);
pause(0.5);

Tail = 50;
Ylim = [0 5];
st_f = 5000;
nd_f = st_f + 80000;
movl.CurrentTime = st_f/movl.FrameRate;
movd.CurrentTime =  st_f/movl.FrameRate;


for i = st_f:1:nd_f
    Gl = readFrame(movl);
    Gd = readFrame(movd);
    subplot(2,2,1);
    imshow(Gl); hold on;
    for f = fs
        plot(xcl{le}(i,f),ycl{le}(i,f),'*','MarkerSize',10,'Color',Cmap(f,:));
        plot(xhl{le}(i,f),yhl{le}(i,f),'.','MarkerSize',15,'Color',Cmap(f,:)*0.5);
    end
   title([Titles_light{le},' ',num2str(i),' Q: ',num2str(qhl{le}(i,f))])
   hold off

   subplot(2,2,3);
   if i>Tail
    plot(0-Tail:0+Tail,Shl{le}(i-Tail:i+Tail,fish_i_for_speed)); 
    ylim(Ylim); xlim([0-Tail 0+Tail]);hold on;
    plot(0-Tail:0+Tail,Scl{le}(i-Tail:i+Tail,fish_i_for_speed),'r');
    hold off;
    box off;
    xlabel('Time from current frame [frames]');
    ylabel('Speed [pix/frame]');
   end
    subplot(2,2,2);
    imshow(Gd);hold on;
    for f = fs
        plot(xcd{de}(i,f),ycd{de}(i,f),'*','MarkerSize',10,'Color',Cmap(f,:));
        plot(xhd{le}(i,f),yhd{le}(i,f),'.','MarkerSize',15,'Color',Cmap(f,:)*0.5);
    end
   title([Titles_dark{de},' ',num2str(i)])
    hold off
       subplot(2,2,4);
   if i>Tail
        plot(0-Tail:0+Tail,Shd{de}(i-Tail:i+Tail,fish_i_for_speed)); 
        ylim(Ylim); xlim([0-Tail 0+Tail]);
        hold on;
        plot(0-Tail:0+Tail,Scd{de}(i-Tail:i+Tail,fish_i_for_speed),'r');
        hold off;
        box off;
        xlabel('Time from current frame  [frames]');
        ylabel('Speed [pix/frame]');
   end
   pause(0.01);
   if MAKEMOVIE && i>Tail
    frame = getframe(gcf);
    writeVideo(obj,frame);
   end
end

if MAKEMOVIE
   close(obj); 
end
%%
% 
% Meds = zeros(2,Nf);
% for j= 1:Nf
%     Meds(1,j) = median(Shl(qsl(:,j),j));
%     Meds(2,j) = median(Shd(qsd(:,j),j));
% end
% 
% figure();
% temp = Shl(:,1); temp(qsl(:,1)) = nan;
% plot(temp);
% 
% figure();
% temp = Shd(:,1); temp(qsd(:,1)) = nan;
% plot(temp);
% 
% histogram(log(Shl));
%% look at the data:

close all;
night_i = 1; % which night to plot


coor = coord{night_i}; % coorl for light night/ coord for dark night
mov = vid_d{night_i}; % vid_l for light night/ vid_d for dark night

fish_i = [1:6]; % number from 1 to N
points_i = [1:6]; % numbers from 1 to n points 
start_f = 10000; % start frame
end_f = start_f+5000; % end frame

plotFishMovie(coor,mov,'fish_i',fish_i,'points_i',points_i,'start_f',start_f,...
    'end_f',end_f,'data_per_point',3,'points_per_fish',Num_dots)

%%
