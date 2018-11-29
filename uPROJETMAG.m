% uPROJETMAG: a simple user interface system for recovering patterns 
% repeating in time and frequency in audio mixtures
% 
% Author: Derry Fitzgerald, Zafar Rafii
% Update: 24/04/2017

function uPROJETMAG

figure(...                                                                  % Create figure window
    'DockControls','off', ...                                               % Interactive figure docking
    'MenuBar','none', ...                                                   % Figure menu bar display
    'Name','uPROJETMAG', ...                                                    % Figure window title
    'NumberTitle','off', ...                                                % Figure window title number
    'Units','normalized', ...                                               % Units of measurement
    'OuterPosition',[0.25,0.25,0.5,0.5]);                                   % Location and size of figure's outer bounds

[A,~,alpha] = imread(fullfile(matlabroot, ...
    'toolbox','matlab','icons','file_open.png'));                           % Read image from graphics file
A = im2double(A);                                                           % Convert image to double precision
A(alpha==0) = NaN;                                                          % Use the AND mask to add transparency
uipushtool( ...                                                             % Create push button on toolbar 
    'CData',A, ...                                                          % Optional image to display on uipushtool
    'ClickedCallback',@open_clickedcallback, ...                            % Callback function that executes when user clicks uipushtool
    'TooltipString','Open an audio file');                                  % Tooltip text

[A,~,alpha] = imread(fullfile(matlabroot, ...
    'toolbox','matlab','icons','file_save.png'));
A = im2double(A);
A(alpha==0) = NaN;
uipushtool( ...
    'CData',A, ...
    'ClickedCallback',@save_clickedcallback, ...
    'TooltipString','Save the audio file');

[A,~,alpha] = imread(fullfile(matlabroot, ...
    'toolbox','matlab','icons','tool_legend.png'));
A = im2double(A);
A(alpha==0) = NaN;
uipushtool( ...
    'CData',A, ...
    'ClickedCallback',@parameters_clickedcallback, ...
    'Separator','on', ...                                                   % Separator line mode
    'TooltipString','Parameters');

addpath('CQT_toolbox_2013')                                                 % Add folder to search path
x = [];                                                                     % Initial input/output
fs = [];
Xcqk = [];
background_or_foreground = 'b';                                             % Initial recovering background (or foreground)
min_time_separation = 1;                                                    % Initial min time separation between repetitions (in seconds)
min_frequency_separation = 1;                                               % Initial min frequency separation between repetitions (in semitones)

    function open_clickedcallback(~,~)
        
        [filename,pathname] = uigetfile( ...                                % Open standard dialog box for retrieving files
            {'*.wav', 'WAVE files (*.wav)'; ...
            '*.mp3', 'MP3 files (*.mp3)'}, ...
            'Select an audio file');
        if isequal(filename,0)                                              % Return if user selects Cancel
            return
        end
        
        file = fullfile(pathname,filename);                                 % Build full file name from parts
        [x,fs] = wavread(file);                                           % Read audio file
        [l,p] = size(x);                                                    % Number of samples and channels
        if p>1
            x=mean(x,2);                                                    %Make file mono for now, do stereo once proof of concept is working
        end
        B = 24;                                                             % Number of bins per octave
        fmin = 27.50;                                                       % Lowest frequency to be analyzed
        fmax = fs/2;                                                        % Highest frequency to be analyzed
                                                     
        V = [];                                                             % Magnitude spectrogram

        Xcqk = cqt(x,B,fs,fmin,fmax);                              % Constant Q Transform (CQT)
        V = abs(Xcqk.c);

        [n,m,~] = size(V);                                                  % Number of frequency channels and time frames
        
         P = db(V);                                                  % Mean spectrogram in decibels
        imagesc(P);                                                         % Scale data and display image object
        set(gca, ...
            'CLim',[min(P(:)),max(P(:))], ...                               % Color limits for objects using colormap
            'XTick',round(m/(l/fs)):round(m/(l/fs)):m, ...                  % X-tick mark locations (every 1 second)
            'XTickLabel',1:round(l/fs), ...                                 % X-tick mark labels (every 1 second)
            'YDir','normal', ...                                            % Direction of increasing values along axis
            'YTick',1:B:n, ...                                              % Y-tick mark locations (every "A" Hz)
            'YTickLabel',fmin*2.^(0:ceil(n/B)-1))                           % Y-tick mark labels (every "A" Hz)
        title(filename, ...                                                 % Add title to current axes
            'Interpreter','none')                                           % Interpretation of text characters
        xlabel('time (s)')                                                  % Label x-axis
        ylabel('log-frequency (Hz)')                                        % Label y-axis
        
        while 1                                                             % Infinite loop
            h = imrect(gca);                                                % Create draggable rectangle
            if isempty(h)                                                   % Return if figure close
                return
            end
            fcn = makeConstrainToRectFcn('imrect', ...
                get(gca,'XLim'),get(gca,'YLim'));                           % Create rectangularly bounded drag constraint function
            setPositionConstraintFcn(h,fcn);                                % Set position constraint function of ROI object
            position = wait(h);                                             % Block MATLAB command line until ROI creation is finished
            if isempty(position)                                            % Return if figure close
                return
            end
            delete(h)                                                       % Remove files or objects
            
            b = waitbar(0,'Please wait...');                                % Open wait bar dialog box
            j = round(position(1));                                         % X-position
            i = round(position(2));                                         % Y-position
            w = round(position(3));                                         % Width
            h = round(position(4));                                         % Height
            R = V(i:i+h-1,j:j+w-1);                                       % Selected rectangle
            C = normxcorr2(R,V);                            % Normalized 2-D cross-correlation
            V = padarray(V,[h-1,w-1],'replicate');                        % Pad array for finding peaks
            
            mpd = [min_frequency_separation*2, ...
                min_time_separation*round(m/(l/fs))];                       % Minimum peak separation
                                                  
                [~,I] = max(C(:));                                          % Linear index of peak
                [I,J] = ind2sub([n+h-1,m+w-1],I);                           % Subscripts from linear index
                C(max(1,I-mpd(1)+1):min(n+h-1,I+mpd(1)+1), ...
                    max(1,J-mpd(2)+1):min(m+w-1,J+mpd(2)+1)) = 0;          % Zero neighborhood around peak
                [~,I] = max(C(:));                                          % Linear index of peak
                [I,J] = ind2sub([n+h-1,m+w-1],I);                           % Subscripts from linear index
                R1 = V(I:I+h-1,J:J+w-1,:);                          % Get most similar rectangle
                %waitbar(k/np,b)                                             % Update wait bar dialog box
                                                                
                [nr,nc]=size(R);
                R=reshape(R,nr*nc,1);                               %Vectorize spectrogram rectangles to perform PROJETMAG
                R1=reshape(R1,nr*nc,1);
           % close(b)                                                        % Close wait bar dialog box
            
            V = V(h:n+h-1,w:m+w-1,:);                                       % Remove pad array
            %Do PROJETMAG on chosen region
            %set up LCR panning (equivalent to assuming repeating pattern
            %is equally loud in both of the regions)
            pan_dir=[0 pi/4 pi/2];
            np=3;
            projmat=[sin(pan_dir); -cos(pan_dir)];

            panmat=[cos(pan_dir); sin(pan_dir)];
            K=abs(panmat'*projmat);
            X=[R R1];
            
            Cc=X*projmat;
            C=abs(Cc);

            %initialise parameters
            %S contains estimated spectrograms of the underlying sources
            S=abs(randn(nr*nc,np))+1;
            lS=repmat(sum(K,2)',nr*nc,1)+eps;
            niter=200;
            %iteratively improve estimates of S
            for k=1:niter
                est=S*K;
                ei=C./(est);
                uS=ei*K';
                S=S.*(uS./lS).^2; %power of 2 to speed up convergence (heuristic)
                waitbar(k/niter,b)
               
            end
            est=S*K;
            close(b);
            %get phase of original excerpt
            fang=cos(angle(Xcqk.c(i:i+h-1,j:j+w-1))) +1i*sin(angle(Xcqk.c(i:i+h-1,j:j+w-1)));
            if strcmp(background_or_foreground,'f')
                reind=1; %foreground of chosen excerpt is in (S(:,1))
            else
                reind=2; %background of chosen excerpt is in (S(:,2))
            end
            ck=C.*((S(:,reind)*K(reind,:)).^2)./((est+eps).^2);
            cka=reshape(ck(:,1),nr,nc);
            %re-apply original phase information
            Xcqk.c(i:i+h-1,j:j+w-1)=cka.*fang;
           
            P = getimage(gca);                                              % Image data from axes
            %P(i:i+h-1,j:j+w-1) = ck;
            P(i:i+h-1,j:j+w-1) = db(cka);                  % Update rectangle in image
            set(get(gca,'Children'),'CData',P)                              % Update image in axes
       end
        
    end

    function save_clickedcallback(~,~)
        
        if isempty(x)                                                       % Return if no input/output
            return
        end
        
        [filename2,pathname] = uiputfile( ...                               % Open standard dialog box for saving files
            {'*.wav', 'WAVE files (*.wav)'; ...
            '*.mp3', 'MP3 files (*.mp3)'}, ...
            'Save the audio file');
        if isequal(filename2,0)                                             % Return if user selects Cancel
            return
        end
     
        x = icqt(Xcqk);

        file = fullfile(pathname,filename2);                                % Build full file name from parts
        wavwrite(x,fs,file);                                               % Write audio file
        
    end

    function parameters_clickedcallback(~,~)
        
        prompt = {'Recovering background (b) or foreground (f):', ...
            'Min time separation between repetitions (in seconds):', ...
            'Min frequency separation between repetitions (in semitones):'};
        dlg_title = 'Parameters';
        num_lines = 1;
        def = {background_or_foreground, ...
            num2str(min_time_separation), ...
            num2str(min_frequency_separation)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);                  % Create and open input dialog box
        if isempty(answer)                                                  % Return if user selects Cancel
            return
        end
        
        background_or_foreground = answer{1};
        min_time_separation = str2double(answer{2});
        min_frequency_separation = str2double(answer{3});
        
    end
end
