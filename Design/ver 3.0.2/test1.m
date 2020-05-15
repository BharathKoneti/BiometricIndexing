function varargout = test1(varargin)

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @test1_OpeningFcn, ...
                       'gui_OutputFcn',  @test1_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end

% --- Executes just before test1 is made visible.
function test1_OpeningFcn(hObject, eventdata, handles, varargin)
    % Choose default command line output for test1
    handles.output = hObject;
    
    global fingerImage;
    global figures;
    global display1;
    global orientArray;
    global actualImage;
    global moments;
    global cluster;
    global currImage;
    global clusterRefName;
    global clustersFoldName;
    global sourceImagesFoldName;
    global currImgDist;
    global cid;
    global cimg;
    global ridgeDisplay;
    global setCompleteSequence;
    global legendreMoments;
    global displayMomentValues;
    global actImage;
    global actImageLoc;
    
    clusterRefName = 'C:\a\clusterVal';
    clustersFoldName = 'C:\a\Clusters';
    sourceImagesFoldName = 'C:\a\Dbase';
    ridgeDisplay = 0;
    setCompleteSequence = 0;
    display1 = 1;
    legendreMoments = 0;
    figures = 0;
    displayMomentValues = 1;
    
   % Update handles structure
    guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = test1_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % Get default command line output from handles structure
    varargout{1} = handles.output;

% --------------------------------------------------------------------
% --------------------------------------------------------------------

function file_Callback(hObject, eventdata, handles)
    function open_Callback(hObject, eventdata, handles)
        global currImage;
        global fingerImage;
        global figures;
        global display1;
        global actImage;
        global setCompleteSequence;
        global actImageLoc;
         
        set(0,'CurrentFigure',test1);
        cab test1
       
        
        [filename, pathname, filterindex] = uigetfile('*', 'Pick a file');
            set(handles.loc,'String',[pathname,'/',filename]); %loc is a reference to hidden text box

         var1= get(handles.loc,'String');
         actImageLoc = var1;
         rgbImage=imread(var1);
         [rowf, colf, dimf ] = size(rgbImage);
         if(dimf == 3)
             rgbImage=rgb2gray(rgbImage);
         end
         
         if(display1 == 1)
            %figure(test1);
            subplot(1,3,2), imhist(rgbImage);
            subplot(1,3,1), imshow(rgbImage),axis tight;
             
            if(setCompleteSequence == 0)
                 set(handles.left1, 'visible',  'on');  
                 set(handles.left1, 'string' , 'Input Image');
                 set(handles.center1, 'visible',  'on');  
                 set(handles.center1, 'string' , 'Histogram');
            end
         end
         if(figures == 1)
            figure('Name','Input Image'), imshow(rgbImage);
            figure('Name','Histogram of Input Image'), imhist(rgbImage);
         end
        
             set(handles.matched1, 'visible',  'on');  
             set(handles.matched1, 'string' , '');
        
         actImage = rgbImage;
         fingerImage = rgbImage;
         guidata(hObject, handles);

% --------------------------------------------------------------------

function histogram_Callback(hObject, eventdata, handles)
    function equalize_Callback(hObject, eventdata, handles)
        global fingerImage;
        global figures;
        global display1;
        global setCompleteSequence;
        
        rgbImage = fingerImage;
        previousImage = fingerImage;
        
        rgbImage=histeq(rgbImage);
        if(display1 == 1)
            	subplot(1,3,2), imshow(rgbImage),axis tight;
                if(setCompleteSequence == 0)
                     set(handles.center1, 'visible',  'on');  
                     set(handles.center1, 'string' , 'Equalized Image');
                end
        end
        if(figures == 1)
            figure('Name','Histogram of Equalized Image'), imhist(rgbImage);
            figure('Name','Equalized Image'), imshow(rgbImage);
        end
        
        
        
        fingerImage=rgbImage;
        guidata(hObject, handles);
        
% --------------------------------------------------------------------

function ridge_Callback(hObject, eventdata, handles) 
    function binarization_Callback(hObject, eventdata, handles)
        global fingerImage;
        global figures;
        global display1;
        global setCompleteSequence;
        
        a = fingerImage;
        previousImage = a;
        BinarizedImage=a(:,:,1)>160;
        if(display1 == 1)
            subplot(1,3,2), imshow(previousImage);
            subplot(1,3,3), imshow(BinarizedImage);
            
            if(setCompleteSequence == 0)
                set(handles.center1, 'visible',  'on');  
                set(handles.center1, 'string' , 'Previous Process');
                set(handles.right1, 'visible',  'on');  
                set(handles.right1, 'string' , 'Binarized Image');
            end
        end
        if(figures == 1)
            figure('name','Binarized Image'), imshow(BinarizedImage);
        end
        
        fingerImage = BinarizedImage;
        guidata(hObject, handles);
    function orientation_Callback(hObject, eventdata, handles)
        global fingerImage;
        global display1;
        global orientArray;
        global ridgeDisplay;
        rgbImage = fingerImage;
        rgbImage = imresize(rgbImage,[335,335],'bilinear');
        gradientMultiplier = 1;        
        orientSmoothMultiplier = 3;
        blockMultiplier = 3;

        [rows,cols] = size(rgbImage);

        blockSize = fix(6*gradientMultiplier);   if ~mod(blockSize,2); blockSize = blockSize+1; end
        filterGauss = fspecial('gaussian', blockSize, gradientMultiplier);
        [fx,fy] = gradient(filterGauss);                       

        gradientX = filter2(fx, rgbImage); 
        gradientY = filter2(fy, rgbImage); 

        gradientXSquare = gradientX.^2;    
        gradientXY = gradientX.*gradientY;
        gradientYSquare = gradientY.^2;

        blockSize = fix(6*blockMultiplier);   if ~mod(blockSize,2); blockSize = blockSize+1; end    
        filterGauss = fspecial('gaussian', blockSize, blockMultiplier);
        gradientXSquare = filter2(filterGauss, gradientXSquare); 
        gradientXY = 2*filter2(filterGauss, gradientXY);
        gradientYSquare = filter2(filterGauss, gradientYSquare);

        divisorVxy = sqrt(gradientXY.^2 + (gradientXSquare - gradientYSquare).^2) + eps;
        sinTheta = gradientXY./divisorVxy; 
        cosTheta = (gradientXSquare-gradientYSquare)./divisorVxy;

        if orientSmoothMultiplier
            blockSize = fix(6*orientSmoothMultiplier);   if ~mod(blockSize,2); blockSize = blockSize+1; end    
            filterGauss = fspecial('gaussian', blockSize, orientSmoothMultiplier);    
            cosTheta = filter2(filterGauss, cosTheta);
            sinTheta = filter2(filterGauss, sinTheta);
        end

        orientim = pi/2 + atan2(sinTheta,cosTheta)/2;
        
        orientArray = orientim;
        
        set(handles.text4, 'visible',  'on');  
        set(handles.text4, 'string' , 'Local Ridge Orientation Estimated');
        
        if(ridgeDisplay == 1)
            plotridgeorient(orientim,10,rgbImage,4); 
        end
        guidata(hObject, handles);
       function plotridgeorient(orient, spacing, im, figno)

        if fix(spacing) ~= spacing
        error('spacing must be an integer');
        end

        %figure('name','orientation image'), imshow(orient);

        [rows, cols] = size(orient);

        lw = 2;      
        len = 0.8*spacing;

        s_orient = orient(spacing:spacing:rows-spacing, ...
                  spacing:spacing:cols-spacing);

        xoff = len/2*cos(s_orient);
        yoff = len/2*sin(s_orient);    

        figure('name','Ridge Image'), imshow(im); hold off
   

        [x,y] = meshgrid(spacing:spacing:cols-spacing, ...
                 spacing:spacing:rows-spacing);

        x = x-xoff;
        y = y-yoff;

        u = xoff*2;
        v = yoff*2;

        quiver(x,y,u,v,0,'.','linewidth',1, 'color','r');

        axis equal, axis ij,  hold on                                 
  
% --------------------------------------------------------------------

function tools_Callback(hObject, eventdata, handles)
    function figure_Callback(hObject, eventdata, handles)
        function figuresOn_Callback(hObject, eventdata, handles)
            global figures;
            figures = 1;
            guidata(hObject, handles);
        function figuresOff_Callback(hObject, eventdata, handles)
            global figures;
            figures = 0;
            guidata(hObject, handles);
    function display_Callback(hObject, eventdata, handles)
        function displayOn_Callback(hObject, eventdata, handles)
            global display1;
            display1 = 1;
            guidata(hObject, handles);
        function displayOff_Callback(hObject, eventdata, handles)
            global display1;
            display1 = 0;
            guidata(hObject, handles);
    function ridgeDisplay_Callback(hObject, eventdata, handles)
        function ridgeOn_Callback(hObject, eventdata, handles)
            global ridgeDisplay;
            ridgeDisplay = 1;
            guidata(hObject, handles);
        function ridgeOff_Callback(hObject, eventdata, handles)
            global ridgeDisplay;
            ridgeDisplay = 0;
            guidata(hObject, handles);
    function momentsD_Callback(hObject, eventdata, handles)
        function momentsOn_Callback(hObject, eventdata, handles)
            global legendreMoments;
            legendreMoments = 1;
            guidata(hObject, handles);
        function momentsOff_Callback(hObject, eventdata, handles)
            global legendreMoments;
            legendreMoments = 0;
            guidata(hObject, handles);
    function displayPolynomialValues_Callback(hObject, eventdata, handles)
        function momentValuesOn_Callback(hObject, eventdata, handles)
            global displayMomentValues;
            displayMomentValues = 1;
            guidata(hObject, handles)
        function momentValuesOff_Callback(hObject, eventdata, handles)
            global displayMomentValues;
            displayMomentValues = 0;
            guidata(hObject, handles)

% --------------------------------------------------------------------

function moments_Callback(hObject, eventdata, handles)
    function legendre_Callback(hObject, eventdata, handles)
        global figures;
        global display1;
        global orientArray;
        global moments;
        global fingerImage;
        global setCompleteSequence;
        global legendreMoments;
        global displayMomentValues;

        
        in = orientArray;     
        rgbImage = fingerImage;
        max = 5;
     

        legMatMap = Legendre_Matrix(in,max);

        in = legMatMap;
        in=fix(in);
        
        if(displayMomentValues == 1)
           in
        end
        if(display1 == 1)
            if(legendreMoments == 1)
                subplot(1,3,2), imshow(rgbImage);
                subplot(1,3,3), imshow(in);
                if(setCompleteSequence == 0)
                    set(handles.center1, 'visible',  'on');  
                    set(handles.center1, 'string' , 'Previous Process');
                    set(handles.right1, 'visible',  'on');  
                    set(handles.right1, 'string' , 'Estimated Moments');
                end
            end
        end
        
        if(figures == 1)
            figure('Name','Legendre Moments'), imshow(in)% hObject    handle to sobel (see GCBO)
        end
        
        set(handles.text4, 'visible',  'on');  
        set(handles.text4, 'string' , 'Legendre Moments Calculated');
       
        
        moments = in; 
        guidata(hObject, handles);
        function coe = Legendre_Matrix(img,max)
            [nx,ny]=size(img);
            N = 128;
            img = double(img);
            X = zeros(nx,ny);
            Y = zeros(nx,ny);

            for x = 0 : nx - 1
                for n = 0 : max
                    if(n == 0)
                        M_LX(n + 1,x + 1) = 1;
                    elseif(n == 1)
                        M_LX(n + 1,x + 1) = (2*x/(N-1)-(nx-1)/(N-1));
                    else 
                        A = (2*n-1)*( 2*x/(N-1)-(nx-1)/(N-1))/n;
                        B = (n-1)/n;
                        M_LX(n + 1,x + 1) = (M_LX(n,x + 1) * A - M_LX(n-1,x + 1) * B);
                    end
                end
            end

            for y = 0 : ny - 1
                for n = 0 : max
                    if(n == 0)
                        M_LY(n + 1,y + 1) = 1;
                    elseif(n == 1)
                        M_LY(n + 1,y + 1) = (2*y/(N-1)-(ny-1)/(N-1));
                    else 
                        A = (2*n-1)*( 2*y/(N-1)-(ny-1)/(N-1))/n;
                        B = (n-1)/n;
                        M_LY(n + 1,y + 1) = (M_LY(n,y + 1) * A - M_LY(n-1,y + 1) * B);
                    end
                end
            end

            for i = 0 : max 
                for j = 0 : max
                        px = M_LX(i + 1, :);
                        py = M_LY(j + 1, :);
                        Q = px' * py;
                        Q = imresize(Q,[335,335],'bilinear');



                        temp = (img .* Q);
                        coe(i + 1,j + 1) = sum(temp(:));
                        coe(i + 1,j + 1) = coe(i + 1,j + 1) * 4 /((N)^2);
                end
            end

% --------------------------------------------------------------------

function clustering_Callback(hObject, eventdata, handles)
    function kmeans_Callback(hObject, eventdata, handles)
        global clustersFoldName;

        dloc = clustersFoldName;
        dloc = strcat(dloc,'\');
        
        a = fullfile(dloc,'cluster1');
        if((exist(a,'dir'))~=7)
            mkdir(fullfile(dloc,'cluster1'));
            mkdir(fullfile(dloc,'cluster2'));
            mkdir(fullfile(dloc,'cluster3'));
            mkdir(fullfile(dloc,'cluster4'));
            mkdir(fullfile(dloc,'cluster5'));
        end
        test1('euclid_Callback',hObject,eventdata,guidata(hObject));
        test1('writeImage_Callback',hObject,eventdata,guidata(hObject));
        function euclid_Callback(hObject, eventdata, handles)
         global clusterRefName;
         global currImgDist;
         global cluster;
         global moments;
         
         in = moments; 
         
         pathname = clusterRefName;
         pathname = strcat(pathname,'\');
         pathnameext = strcat(pathname,'*.m');
         fnames = dir(pathnameext);
         numfids = length(fnames);
         small = 100000;
         clusterIndex = 0;
         for k = 1:numfids
             floc = [pathname,fnames(k).name];
             fid = fopen(floc);
                tline = fgetl(fid);
             fclose(fid);
             ref_array2 = str2num(tline);
             [m1,n1]=size(ref_array2);
             [m2,n2]=size(in);
             D=zeros(m1,m2);
           
             for i=1:m1
               for j=1:m2
                    D(i,j)=sqrt((ref_array2(i,1)-in(j,1))^2+(ref_array2(i,2)-in(j,2))^2);
                end
             end
            
             D=fix(D);

             avgVal = sum(D)/6;
             if(avgVal < small)
                 small = avgVal;
                 clusterIndex = k;
             end
         end
         
         currImgDist = small;
         
         cluster = clusterIndex;
         guidata(hObject, handles);
    function locateImageCluster_Callback(hObject, eventdata, handles)
        global clustersFoldName;
        global cluster;
        global currImgDist;
        global cid;

        test1('euclid_Callback',hObject,eventdata,guidata(hObject));
        
        clusterFolder = cluster;
        
        cid = currImgDist;
        
        dloc =  clustersFoldName;
        dloc = strcat(dloc,'\');
        
        f1 = strcat(dloc,'Cluster1\');
        f2 = strcat(dloc,'Cluster2\');
        f3 = strcat(dloc,'Cluster3\');
        f4 = strcat(dloc,'Cluster4\');
        f5 = strcat(dloc,'Cluster5\');
        global sourceImagesFoldName;
        
        switch(clusterFolder)
            case 1 
                sourceImagesFoldName = f1;
            case 2 
                sourceImagesFoldName = f2;
            case 3 
                sourceImagesFoldName = f3;
            case 4
                sourceImagesFoldName = f4;
            case 5
                sourceImagesFoldName = f5;
        end
        
        global setCompleteSequence; 
        setCompleteSequence = 0;
        guidata(hObject, handles);
        
        test1('orderSequenceNonStop_Callback',hObject,eventdata,guidata(hObject));
                
             
% --------------------------------------------------------------------

function executeInOrder_Callback(hObject, eventdata, handles)
    function executeOrderSequence_Callback(hObject, eventdata, handles)
         global setCompleteSequence;
         setCompleteSequence = 1;
         test1('orderSequenceNonStop_Callback',hObject,eventdata,guidata(hObject));
         test1('displayOn_Callback',hObject,eventdata,guidata(hObject))
        function orderSequenceNonStop_Callback(hObject, eventdata, handles)
            global sourceImagesFoldName;
            global clustersFoldName;
            global setCompleteSequence;
            global actualImage;
            global fingerImage;
            global currImage;
            global currImgDist;
            global cid;
            global figures;
            
             
            test1('momentValuesOff_Callback',hObject,eventdata,guidata(hObject));
            test1('displayOff_Callback',hObject,eventdata,guidata(hObject));
            test1('figuresOff_Callback',hObject,eventdata,guidata(hObject));
            imgPathNameRaw = sourceImagesFoldName;
            imgpathname = strcat(imgPathNameRaw,'\');
            imgdirname = strcat(imgpathname,'*.TIF');
            imgnames = dir(imgdirname);
            imgfids = length(imgnames);

           
            dloc =  clustersFoldName;
            
            dloc = strcat(dloc,'\');

            a = fullfile(dloc,'cluster1');
            if((exist(a,'dir'))~=7)
                mkdir(fullfile(dloc,'cluster1'));
                mkdir(fullfile(dloc,'cluster2'));
                mkdir(fullfile(dloc,'cluster3'));
                mkdir(fullfile(dloc,'cluster4'));
                mkdir(fullfile(dloc,'cluster5'));
            end
            
            count = 0;
            for z = 1:imgfids
                str1 = strcat(imgpathname,imgnames(z).name);
                rgbImage=imread(str1);
                 [rowf, colf, dimf ] = size(rgbImage);
                 if(dimf == 3)
                     rgbImage=rgb2gray(rgbImage);
                 end
                 %imgnames(z).name
                 cimg = rgbImage;
                 actualImage = rgbImage;
                 fingerImage = rgbImage;

                 currImage = str1;
                 guidata(hObject, handles);

                
                
                test1('equalize_Callback',hObject,eventdata,guidata(hObject));
                test1('binarization_Callback',hObject,eventdata,guidata(hObject));
                test1('orientation_Callback',hObject,eventdata,guidata(hObject));
                test1('legendre_Callback',hObject,eventdata,guidata(hObject));
                test1('euclid_Callback',hObject,eventdata,guidata(hObject));
                if(setCompleteSequence == 0)
                    vcid = currImgDist;
                    acid = cid;
                    
                    displayCl = strcat('Images are being recognized from the cluster : ',sourceImagesFoldName);
                    
                    set(handles.text4, 'visible',  'on');  
                    set(handles.text4, 'string' , displayCl);
                    
                    if(abs(acid-vcid)==0)
                         figure('Name','Matched Image'), imshow(actualImage);
                         count = count + 1;
                    else if( (abs(acid-vcid)<=15) )
                         figure('Name','Nearly Matched Image'), imshow(actualImage);
                        count = count + 1;
                        end
                    end
                    
                end
                if(setCompleteSequence == 1)
                    test1('writeImage_Callback',hObject,eventdata,guidata(hObject));
                end
            end
            if(setCompleteSequence == 0)
                if(count == 0)
                    set(handles.matched1, 'visible',  'on');  
                    set(handles.matched1, 'string' , 'No Matches Found');
                end
            end
            test1('momentValuesOn_Callback',hObject,eventdata,guidata(hObject));
        function writeImage_Callback(hObject, eventdata, handles)
            global currImage;
            global clustersFoldName;
            global cluster;
            global actImage;
            global actualImage;
            global setCompleteSequence;
            global actImageLoc;
            
            c1 = cluster;
            if(setCompleteSequence == 1)
                image = actualImage;
                [imgPath, imgName, imgExt] = fileparts(currImage);
            else
                image = actImage;
                [imgPath, imgName, imgExt] = fileparts(actImageLoc);
            end
            
            
            dloc =  clustersFoldName;
            dloc = strcat(dloc,'\');
            
            dispCluster = strcat('Fingerprint sorted to Cluster',num2str(c1),' folder under Location:   ',dloc); 
            set(handles.text4, 'visible',  'on');  
            set(handles.text4, 'string' , dispCluster);
           
            
            f1 = strcat(dloc,'Cluster1\');
            f2 = strcat(dloc,'Cluster2\');
            f3 = strcat(dloc,'Cluster3\');
            f4 = strcat(dloc,'Cluster4\');
            f5 = strcat(dloc,'Cluster5\');
            
            s = strcat(imgName,imgExt);

            switch(c1)
                case 1
                    imwrite(image, fullfile(f1, s));
                case 2
                    imwrite(image, fullfile(f2, s));
                case 3
                    imwrite(image, fullfile(f3, s));
                case 4
                    imwrite(image, fullfile(f4, s));
                case 5
                    imwrite(image, fullfile(f5, s));
            end
            
% --------------------------------------------------------------------

function environmentVariables_Callback(hObject, eventdata, handles)
    function clusterVariables_Callback(hObject, eventdata, handles)
        function clusterReferenceValues_Callback(hObject, eventdata, handles)
            global clusterRefName;
            clusterRefName = uigetdir;
            guidata(hObject, handles);           
        function clustersFolders_Callback(hObject, eventdata, handles)
            global clustersFoldName;
            clustersFoldName = uigetdir;
            guidata(hObject, handles);
        function sourceImageDatabase_Callback(hObject, eventdata, handles)
           global sourceImagesFoldName;
           sourceImagesFoldName = uigetdir;
           guidata(hObject, handles);

% --------------------------------------------------------------------
