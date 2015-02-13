function StimulusPresentationRT(bugMe)

%% Imports
import Devices.NI.DAQmx.*

%% Initiations and instatiations of persistent variables
counter=0;
audio_freq = 20000;
wrong = 0;
afterLick = 0;
persistent hCtr data
persistent w rect frameRate
persistent dotCentre dotCentrePolar dotDirection
persistent elapsedTime lickTime
persistent player
persistent s1 s2

close all;

%% Interrupt section
if nargin < 1
    killMe = false;
elseif bugMe
    killMe = bugMe;
end

if killMe && ~isempty(hCtr)
    hCtr.stop();
    return;
elseif ~isempty(hCtr)
    delete(hCtr);
end

%% Opens Psychtoolbox window
% Open a double buffered fullscreen window and select a gray background
% color:
screenNumber = 2; % Stimulus screen
if isempty(w)
    [w, rect]  = Screen('OpenWindow', screenNumber, 0,[], 8, 2);
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    frameRate=Screen('FrameRate',screenNumber);
    
end

%% NIDAQ Configuration
%Create AI Task
deviceName = 'Dev1';
hCtr = Task('myTask');
hCtr.createAIVoltageChan(deviceName,0,'myTriggerAI0');
hCtr.cfgDigEdgeStartTrig('PFI0');
hCtr.cfgSampClkTiming(10000,'DAQmx_Val_FiniteSamps',2);
hCtr.registerDoneEvent(@FrameStartCallback); %Both DAQmx_Val_SampleClock and DAQmx_Val_SampleCompleteEvent yield same results
pause(1); %Give the registration a chance to 'settle' before starting Task. (otherwise, first event is skipped)
hCtr.start();


%% Trigger received
    function FrameStartCallback(~,~)
        %disp('Trigger received')
        switch counter
            case 0
                counter = 1;
                PrepareStimulus();
                stop(hCtr);
                start(hCtr);
                
            case 1
                counter = 2;
                stop(hCtr);
                start(hCtr);
                PresentStimulus();
                counter=0;
                
            case 2
                counter = 3;
                LickDetected();
                stop(hCtr);
                start(hCtr);
                
            case 3
                counter = 4;
                TimerOut();
                stop(hCtr);
                start(hCtr);
            case 4
                AfterLick();
                stop(hCtr);
                start(hCtr);
        end
        
    end

    function PrepareStimulus()
        data = load('c:/ratter/next_trial');
        afterLick = 0;
        
        s1 = RandStream('mlfg6331_64','seed',data.rand_seed);
        s2 = RandStream('mlfg6331_64','seed',0);
        
        % Dots speed
        dotSpeedPerFrame = data.dot_speed_px/frameRate;
        
        % Initiates the dots' centres
        dotCentrePolar = [rand(s1,1,data.dot_number)*2*pi-pi;...
            sqrt(rand(s1,1,data.dot_number))*data.radius_px];
        
        nDotLimit = round(data.dot_number*(1-data.dot_coherence)/8);
        sDotLimit = round(data.dot_number*(1-data.dot_coherence)/4);
        eDotLimit = round(data.dot_number*(1-data.dot_coherence)*3/8);
        wDotLimit = round(data.dot_number*(1-data.dot_coherence)/2);
        neDotLimit = round(data.dot_number*(1-data.dot_coherence)*5/8);
        seDotLimit = round(data.dot_number*(1-data.dot_coherence)*3/4);
        swDotLimit = round(data.dot_number*(1-data.dot_coherence)*7/8);
        nwDotLimit = round(data.dot_number*(1-data.dot_coherence));
        
        dotDirection = zeros(2,data.dot_number);
        
        % Change the way direction and orientation is assessed
        % Just two orientations and rotate the whole frame by a specific
        % angle
        dotDirection(:,1:nDotLimit) = repmat([0; -1]*dotSpeedPerFrame,1,nDotLimit);
        dotDirection(:,nDotLimit+1:sDotLimit) = repmat([0; 1]*dotSpeedPerFrame,1,sDotLimit-nDotLimit);
        dotDirection(:,sDotLimit+1:eDotLimit) = repmat([1; 0]*dotSpeedPerFrame,1,eDotLimit -sDotLimit);
        dotDirection(:,eDotLimit+1:wDotLimit) = repmat([-1; 0]*dotSpeedPerFrame,1,wDotLimit-eDotLimit);
        dotDirection(:,wDotLimit+1:neDotLimit) = repmat([sqrt(2)/2; -sqrt(2)/2]*dotSpeedPerFrame,1,neDotLimit-wDotLimit);
        dotDirection(:,neDotLimit+1:seDotLimit) = repmat([sqrt(2)/2; sqrt(2)/2]*dotSpeedPerFrame,1,seDotLimit-neDotLimit);
        dotDirection(:,seDotLimit+1:swDotLimit) = repmat([-sqrt(2)/2; sqrt(2)/2]*dotSpeedPerFrame,1,swDotLimit-seDotLimit);
        dotDirection(:,swDotLimit+1:nwDotLimit) = repmat([-sqrt(2)/2; -sqrt(2)/2]*dotSpeedPerFrame,1,nwDotLimit-swDotLimit);
        
        if data.dot_coherence>0
            coherentDirection = [sin(data.stim_dir*pi/180); -cos(data.stim_dir*pi/180)];
            
            dotDirection(:,nwDotLimit+1:end) =  repmat(coherentDirection*dotSpeedPerFrame,1,data.dot_number-nwDotLimit);
        end
        
        disp('stimulus prepared');
        
        
    end

    function PresentStimulus()
        
        nFramesStim = round(data.stim_length*frameRate);
        nFrames = round((data.resp_delay+data.resp_window+data.after_lick+0.3)*frameRate);
        
        priorityLevel=MaxPriority(w);
        Priority(priorityLevel);
        
        % Maybe include permutation of dotDirection so that the
        % order of position reset is randomized over the different
        % directions
        dotRem = rem(1:data.dot_number,data.dot_lifetime);
        
        Screen('FillRect',w, uint8(data.background_level),rect);
        
        if data.sound_sides
            if data.left
                player = audioplayer([zeros(1, data.pre_stimulus*audio_freq),sin((1:audio_freq*(data.stim_length))/audio_freq*2*pi*(data.left_frequency)*1000)]*(data.sound_sides_volume)/2000,audio_freq);
            else
                player = audioplayer([zeros(1, data.pre_stimulus*audio_freq),sin((1:audio_freq*(data.stim_length))/audio_freq*2*pi*(data.right_frequency)*1000)]*(data.sound_sides_volume)/2000,audio_freq);
            end
            play(player);
        end
        
        disp('stimulus started');
        
        for i=1:nFrames
            if afterLick
                afterLick = 0;
                if not(wrong)
                    disp('correct lick')
                    if data.trial_sounds
                        player = audioplayer(sin((1:audio_freq*(data.reward_sound_length))/audio_freq*2*pi*1000)*(data.sound_volume)/1000,audio_freq);
                        play(player);
                    end
                else
                    disp('error lick')
                    if data.trial_sounds
                        player = audioplayer((rand(s2,audio_freq*(data.error_sound_length),1) - 0.5)*(data.sound_volume)/1000,audio_freq);
                        play(player);
                    end
                end
                
                if data.stop_stimulus && ...
                        (not(wrong) || (wrong && not(data.visual_error)))
                    Screen('FillRect',w, uint8(data.inter_stim_level),rect);
                    Screen('Flip', w);
                    Priority(0);
                    disp('stimulus stoped');
                    return;
                elseif  wrong && data.visual_error
                    nFramesError = round(data.error_length*frameRate);
                    for j=1:nFramesError
                        if mod(j,data.error_lifetime)/data.error_lifetime<=0.5
                            Screen('FillRect',w, uint8(data.background_level),rect);
                        else
                            Screen('FillRect',w, uint8(data.stim_level)*0.2,rect);
                        end
                        Screen('Flip', w);
                    end
                    Screen('FillRect',w, uint8(data.inter_stim_level),rect);
                    Screen('Flip', w);
                    Priority(0);
                    disp('error stimulus')
                    return;
                end
            end
            
            if i<nFramesStim
                % Converts dots positions to cartisian coordinates
                dotCentre = [cos(dotCentrePolar(1,:)).*dotCentrePolar(2,:); sin(dotCentrePolar(1,:)).*dotCentrePolar(2,:)];
                
                Screen('DrawDots', w, dotCentre, data.dot_size_px, uint8(data.stim_level), data.centre_px,1);
                Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                Screen('Flip', w);
                
                % Moves dots
                dotCentre = dotCentre + dotDirection;
                
                % Converts updated dots positions to polar coordinates
                dotCentrePolar = [atan2(dotCentre(2,:),dotCentre(1,:)); sqrt(dotCentre(1,:).^2 + dotCentre(2,:).^2)];
                
                % Redraws the position of dead dots
                frameRem = rem(i,data.dot_lifetime);
                dotCentrePolar(:,dotRem==frameRem) = ...
                    [rand(s1,1,sum(dotRem==frameRem))*2*pi-pi;...
                    sqrt(rand(s1,1,sum(dotRem==frameRem)))*data.radius_px];
                
                % Detects whether a border was reached
                dotCentrePolar(1, dotCentrePolar(2,:) > data.radius_px) = dotCentrePolar(1, dotCentrePolar(2,:) > data.radius_px)+pi;
            elseif i == nFramesStim
                disp('stimulus presented');
                Screen('FillRect',w, uint8(data.inter_stim_level),rect);
                Screen('Flip', w);
            else
                Screen('Flip', w);
            end
            
        end
        Priority(0);
    end

    function LickDetected()
        disp('lick');
        lickTime = tic;
    end

    function TimerOut()
        elapsedTime = toc(lickTime);
        
        if elapsedTime < 0.15
            wrong = 0;
        else
            wrong = 1;
        end
        
        if data.after_lick==0
            afterLick = 1;
        end
    end

    function AfterLick()
        afterLick = 1;
    end

end
