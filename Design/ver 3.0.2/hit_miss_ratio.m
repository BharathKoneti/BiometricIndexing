function varargout = hit_miss_ratio(varargin)
% HIT_MISS_RATIO MATLAB code for hit_miss_ratio.fig
%      HIT_MISS_RATIO, by itself, creates a new HIT_MISS_RATIO or raises the existing
%      singleton*.
%
%      H = HIT_MISS_RATIO returns the handle to a new HIT_MISS_RATIO or the handle to
%      the existing singleton*.
%
%      HIT_MISS_RATIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HIT_MISS_RATIO.M with the given input arguments.
%
%      HIT_MISS_RATIO('Property','Value',...) creates a new HIT_MISS_RATIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hit_miss_ratio_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hit_miss_ratio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hit_miss_ratio

% Last Modified by GUIDE v2.5 07-Apr-2013 19:28:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hit_miss_ratio_OpeningFcn, ...
                   'gui_OutputFcn',  @hit_miss_ratio_OutputFcn, ...
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

function hit_miss_ratio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hit_miss_ratio (see VARARGIN)

% Choose default command line output for hit_miss_ratio
handles.output = hObject;
global clusterRefName;

clusterRefName = 'C:\a\clusterVal';
% Update handles structure
guidata(hObject, handles);

function varargout = hit_miss_ratio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
% --------------------------------------------------------------------

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
        
function sdb_Callback(hObject, eventdata, handles)
    dir1=uigetdir();
    set(handles.sdt,'String',dir1);
	
    imgpathname = strcat(dir1,'\');
    imgdirname = strcat(imgpathname,'*.TIF');
    imgnames = dir(imgdirname);
    imgfids = length(imgnames);
    set(handles.text5,'String',imgfids);
    
    guidata(hObject, handles);

function dcb_Callback(hObject, eventdata, handles)
    dloc=uigetdir();
    set(handles.dct,'String',dloc);
    
    dloc = strcat(dloc,'\');
    
    imgfids = 0;

    f1 = strcat(dloc,'Cluster1\');
    imgdirname = strcat(f1,'*.TIF');
    imgnames = dir(imgdirname);
    imgfids = imgfids + length(imgnames);
    
    f2 = strcat(dloc,'Cluster2\');
    imgdirname = strcat(f2,'*.TIF');
    imgnames = dir(imgdirname);
    imgfids = imgfids + length(imgnames);
    
    f3 = strcat(dloc,'Cluster3\');
    imgdirname = strcat(f3,'*.TIF');
    imgnames = dir(imgdirname);
    imgfids = imgfids + length(imgnames);
    
    f4 = strcat(dloc,'Cluster4\');
    imgdirname = strcat(f4,'*.TIF');
    imgnames = dir(imgdirname);
    imgfids = imgfids + length(imgnames);
    
    f5 = strcat(dloc,'Cluster5\');
    imgdirname = strcat(f5,'*.TIF');
    imgnames = dir(imgdirname);
    imgfids = imgfids + length(imgnames);
   
    
    set(handles.text6,'String',imgfids);
    
	guidata(hObject, handles);
    
function calculateMatchPercentage_Callback(hObject, eventdata, handles)
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
function orderSequence_Callback(hObject, eventdata, handles)
            global sourceImagesFoldName;
            global clustersFoldName;
            global actualImage;
            global fingerImage;
            global currImage;
            global currImgDist;
            global cid;
            
            sourceImagesFoldName = get(handles.sdt,'String');
            imgPathNameRaw = sourceImagesFoldName;
            imgpathname = strcat(imgPathNameRaw,'\');
            imgdirname = strcat(imgpathname,'*.TIF');
            imgnames = dir(imgdirname);
            imgfids = length(imgnames);

           
            dloc =  get(handles.dct,'String');

            clustersFoldName =  dloc;        
            count = 0;
            for z = 1:imgfids
                str1 = strcat(imgpathname,imgnames(z).name);
                rgbImage=imread(str1);
                 [rowf, colf, dimf ] = size(rgbImage);
                 if(dimf == 3)
                     rgbImage=rgb2gray(rgbImage);
                 end
                 imgnames(z).name
                 actualImage = rgbImage;
                 fingerImage = rgbImage;

                 currImage = str1;
                 guidata(hObject, handles);

                
                
                test1('equalize_Callback',hObject,eventdata,guidata(hObject));
                test1('binarization_Callback',hObject,eventdata,guidata(hObject));
                test1('orientation_Callback',hObject,eventdata,guidata(hObject));
                test1('legendre_Callback',hObject,eventdata,guidata(hObject));
                test1('euclid_Callback',hObject,eventdata,guidata(hObject));

                vcid = currImgDist;
                acid = cid;
                    
                if(abs(acid-vcid)==0)
                    count = count + 1;
                    flag = 1;
                end
            end
