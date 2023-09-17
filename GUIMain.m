function varargout = GUIMain(varargin)
%GUIMAIN M-file for GUIMain.fig
%      GUIMAIN, by itself, creates a new GUIMAIN or raises the existing
%      singleton*.
%
%      H = GUIMAIN returns the handle to a new GUIMAIN or the handle to
%      the existing singleton*.
%
%      GUIMAIN('Property','Value',...) creates a new GUIMAIN using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GUIMain_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUIMAIN('CALLBACK') and GUIMAIN('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUIMAIN.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES 

% Edit the above text to modify the response to help GUIMain

% Last Modified by GUIDE v2.5 06-Aug-2018 15:42:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIMain_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIMain_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUIMain is made visible.
function GUIMain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for GUIMain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIMain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIMain_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Initialize.
function Initialize_Callback(hObject, eventdata, handles)
% hObject    handle to Initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load('Initializer.mat');
set(findobj('Tag','rowem'),'String',rowem); set(findobj('Tag','Eem0'),'String',Eem0);
set(findobj('Tag','ytaem'),'String',ytaem); set(findobj('Tag','sigmaem'),'String',sigmaem);
set(findobj('Tag','rowvm'),'String',rowvm); set(findobj('Tag','Evm0'),'String',Evm0);
set(findobj('Tag','ytavm'),'String',ytavm); set(findobj('Tag','sigmavm'),'String',sigmavm);
set(findobj('Tag','row1'),'String',row1); set(findobj('Tag','c1'),'String',c1);
set(findobj('Tag','row2'),'String',row2);  set(findobj('Tag','c2'),'String',c2);

set(findobj('Tag','dem'),'String',dem); set(findobj('Tag','dvm'),'String',dvm);
set(findobj('Tag','a_inner'),'String',a_inner); set(findobj('Tag','b_middle'),'String',b_middle);
set(findobj('Tag','c_outer'),'String',c_outer);
set(findobj('Tag','fa'),'String',fa); set(findobj('Tag','df'),'String',df);
set(findobj('Tag','fb'),'String',fb);
set(findobj('Tag','cpa'),'String',cpa); set(findobj('Tag','dcp'),'String',dcp);
set(findobj('Tag','cpb'),'String',cpb);
set(findobj('Tag','kia'),'String',kia); set(findobj('Tag','dki'),'String',dki);
set(findobj('Tag','kib'),'String',kib);
set(findobj('Tag','kur'),'String',kur); set(findobj('Tag','err'),'String',err);

switch szFun
    case 'PlanarPlate'
        val = 1; set(findobj('Tag','szFun'),'value',val);
    case 'CylindricalShell'
        val = 2; set(findobj('Tag','szFun'),'value',val);
end
switch szBC
    case 'Va-So-Va'
        val = 1; set(findobj('Tag','szBC'),'value',val);
    case 'Va-So-Fl'
        val = 2; set(findobj('Tag','szBC'),'value',val);
    case 'Fl-So-Fl'
        val = 3; set(findobj('Tag','szBC'),'value',val);
    case 'Va-So-So-Va'
        val = 4; set(findobj('Tag','szBC'),'value',val);
    case 'Va-So-So-Fl'
        val = 5; set(findobj('Tag','szBC'),'value',val);
    case 'Fl-So-So-Fl'
        val = 6; set(findobj('Tag','szBC'),'value',val);
end
switch szMode
    case 'S/L'
        val = 1; set(findobj('Tag','szMode'),'value',val);
    case 'A/T'
        val = 2; set(findobj('Tag','szMode'),'value',val);
    case 'F'
        val = 3; set(findobj('Tag','szMode'),'value',val);
end
set(findobj('Tag','nMode_cs'),'String',nMode_cs);

% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
row1 = str2double(get(findobj('Tag','row1'),'String')); c1 = str2double(get(findobj('Tag','c1'),'String'));
rowem = str2double(get(findobj('Tag','rowem'),'String')); Eem0 = str2double(get(findobj('Tag','Eem0'),'String'));
ytaem = str2double(get(findobj('Tag','ytaem'),'String')); sigmaem = str2double(get(findobj('Tag','sigmaem'),'String'));
rowvm = str2double(get(findobj('Tag','rowvm'),'String')); Evm0 = str2double(get(findobj('Tag','Evm0'),'String'));
ytavm = str2double(get(findobj('Tag','ytavm'),'String')); sigmavm = str2double(get(findobj('Tag','sigmavm'),'String'));
row2 = str2double(get(findobj('Tag','row2'),'String')); c2 = str2double(get(findobj('Tag','c2'),'String'));

dem = str2double(get(findobj('Tag','dem'),'String')); dvm = str2double(get(findobj('Tag','dvm'),'String'));
a_inner = str2double(get(findobj('Tag','a_inner'),'String')); b_middle = str2double(get(findobj('Tag','b_middle'),'String'));
c_outer = str2double(get(findobj('Tag','c_outer'),'String'));
fa = str2double(get(findobj('Tag','fa'),'String')); df = str2double(get(findobj('Tag','df'),'String'));
fb = str2double(get(findobj('Tag','fb'),'String')); cpa = str2double(get(findobj('Tag','cpa'),'String'));
dcp = str2double(get(findobj('Tag','dcp'),'String')); cpb = str2double(get(findobj('Tag','cpb'),'String'));
kia = str2double(get(findobj('Tag','kia'),'String')); dki = str2double(get(findobj('Tag','dki'),'String'));
kib = str2double(get(findobj('Tag','kib'),'String')); 

str = get(findobj('Tag','szFun'),'String'); val = get(findobj('Tag','szFun'),'value'); szFun = str{val};
str = get(findobj('Tag','szBC'),'String'); val = get(findobj('Tag','szBC'),'value'); szBC = str{val};
str = get(findobj('Tag','szMode'),'String'); val = get(findobj('Tag','szMode'),'value'); szMode = str{val};
nMode_cs = str2double(get(findobj('Tag','nMode_cs'),'String'));
kur = str2double(get(findobj('Tag','kur'),'String')); err = str2double(get(findobj('Tag','err'),'String'));

fid= fopen('Initializer.mat');
if fid==-1, readme = 'Parameters File'; save('Initializer.mat','readme'); end
readme = 'Parameters File'; save('Initializer.mat','readme','-append');
save('Initializer.mat','row1','c1','rowem','Eem0','ytaem','sigmaem','rowvm','Evm0','ytavm','sigmavm','row2','c2','-append');
save('Initializer.mat','dem','dvm','a_inner','b_middle','c_outer','-append');
save('Initializer.mat','fa','df','fb','cpa','dcp','cpb','kia','dki','kib','-append');
save('Initializer.mat','szFun','szBC','szMode','nMode_cs','kur','err','-append');
if fid ~= -1; fclose(fid); end

% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc; tic; diary diary.log;
disp('Initializing...');
szCheck = questdlg('To check roots manually?', 'Choice', 'Yes', 'No', 'Yes'); pause(0.1);
%% %参数设定【Parameters Setting】
global szFun szMethod szBC szMode nMode_cs szChkCho
szChkCho = 'Yes';
szMethod = 'Mixed';
%------------------------------------------------------------------------------
rowem = str2double(get(findobj('tag','rowem'),'string'));
Eem0 = str2double(get(findobj('tag','Eem0'),'string'));
ytaem = str2double(get(findobj('tag','ytaem'),'string')); 
sigmaem = str2double(get(findobj('tag','sigmaem'),'string'));
rowvm = str2double(get(findobj('tag','rowvm'),'string')); 
Evm0 = str2double(get(findobj('tag','Evm0'),'string'));
ytavm = str2double(get(findobj('tag','ytavm'),'string')); 
sigmavm = str2double(get(findobj('tag','sigmavm'),'string'));
row1 = str2double(get(findobj('tag','row1'),'string')); 
c1 = str2double(get(findobj('tag','c1'),'string'));
row2 = str2double(get(findobj('tag','row2'),'string')); 
c2 = str2double(get(findobj('tag','c2'),'string'));

dem = str2double(get(findobj('tag','dem'),'string'));
dvm = str2double(get(findobj('tag','dvm'),'string'));
a_inner = str2double(get(findobj('tag','a_inner'),'string'));
b_middle = str2double(get(findobj('tag','b_middle'),'string'));
c_outer = str2double(get(findobj('tag','c_outer'),'string'));

h = findobj('tag','szFun'); str = get(h,'string'); val = get(h,'value'); szFun = str{val};
h = findobj('tag','szBC'); str = get(h,'string'); val = get(h,'value'); szBC = str{val};
nMode_cs = str2double(get(findobj('tag','ncs'),'string'));
kur = str2double(get(findobj('tag','kur'),'string'));
h = findobj('tag','szMode'); str = get(h,'string'); val = get(h,'value'); szMode = str{val};

fa = str2double(get(findobj('tag','fa'),'string'));  
df = str2double(get(findobj('tag','df'),'string')); 
fb = str2double(get(findobj('tag','fb'),'string')); 
cpa = str2double(get(findobj('tag','cpa'),'string'));
dcp = str2double(get(findobj('tag','dcp'),'string')); 
cpb = str2double(get(findobj('tag','cpb'),'string'));
kia = str2double(get(findobj('tag','kia'),'string'));
dki = str2double(get(findobj('tag','dki'),'string'));
kib = str2double(get(findobj('tag','kib'),'string')); 
err = str2double(get(findobj('tag','err'),'string'));
%%
% %参数设定
F = fa:df:fb;
Evm = Evm0*(1-1i*ytavm);
lamdavm = Evm*sigmavm/((1+sigmavm)*(1-2*sigmavm)); miuvm = Evm/(2*(1+sigmavm));
clvm = sqrt((lamdavm+2*miuvm)/rowvm); ctvm = sqrt(miuvm/rowvm);
crvm = RayleighWave(clvm, ctvm);
Eem = Eem0*(1-1i*ytaem);
lamdaem = Eem*sigmaem/((1+sigmaem)*(1-2*sigmaem)); miuem = Eem/(2*(1+sigmaem));
clem = sqrt((lamdaem+2*miuem)/rowem); ctem = sqrt(miuem/rowem);
crem = RayleighWave(clem, ctem); %瑞利波速
%参数保存【Parameters Save】
save('data.mat', 'row1', 'c1', 'row2', 'c2');
save('data.mat', 'rowvm', 'Evm0', 'ytavm', 'sigmavm', 'rowem', 'Eem0', 'ytaem', 'sigmaem', '-append');
save('data.mat', 'lamdavm', 'miuvm', 'clvm', 'ctvm', 'lamdaem', 'miuem', 'clem', 'ctem', '-append');
save('data.mat', 'dem', 'dvm', 'a_inner', 'b_middle', 'c_outer', '-append');
save('data.mat', 'F', 'szFun', 'szBC', 'szMode', '-append');
%% %频散方程求根【Roots Find of Dispersion Equation】
%材料参数【Material Parameters】
MPM = [row1, c1, rowvm, lamdavm, miuvm, dvm, rowem, lamdaem, miuem, dem, row2, c2,...
    a_inner, b_middle, c_outer];
%奇异相速度【Singular Phase Velocities】
SVM = [clem, ctem, crem, clvm, ctvm, crvm];
%搜索参数【Search Parameters】
SPM = [cpa, dcp, cpb, kia, dki, kib];
save('data.mat', 'MPM', 'SVM', 'SPM', '-append');
%%
nL = 1E1;
cp_ps0 = zeros(length(F), nL); cp_ps0(:,:) = NaN; ki_ps0 = cp_ps0;
cp_ps1 = cp_ps0; ki_ps1 = cp_ps0; cp_ps2 = cp_ps0;
ki_ps2 = cp_ps0; cp_ps3 = cp_ps0; ki_ps3 = cp_ps0;
MSW = zeros(length(F), 6); MSW(:, :) = NaN;
n = 1;
disp(['toc: ', num2str(toc), 's']);
disp('*******************************************************');
for f = F
    diary on;
    disp(['f = ', num2str(f), 'Hz']);
    w = 2*pi*f;
    %% %全局搜索【Global Search】
    [cp1, ki1]= GlobalDomainSearch(f, MPM, SPM, kur, err);
    %% %局部搜索【Local Search】
    [cp2, ki2, SW]= LocalDomainSearch(f, MPM, SPM, SVM, kur, err);
    MSW(n, 1:length(SW)) = w./real(SW)+1i*imag(SW); %奇异波数：实部为波速，虚部为衰减系数
    %【Singular Wavenumbers(SW)：real->velocity, imag.->attenuation】
    %% %根值合并【Roots Merge】
    [cp_0, ki_0, cp_1, ki_1, cp_2, ki_2] = RootsMerge(f, SW, cp1, ki1, cp2, ki2);
    cp_ps0(n, 1:length(cp_0)) = cp_0; ki_ps0(n, 1:1:length(ki_0)) = ki_0;
    cp_ps1(n, 1:1:length(cp_1)) = cp_1; ki_ps1(n, 1:1:length(ki_1)) = ki_1;
    cp_ps2(n, 1:1:length(cp_2)) = cp_2; ki_ps2(n, 1:1:length(ki_2)) = ki_2;
    save('data.mat', 'cp_ps0', 'ki_ps0', 'cp_ps1', 'ki_ps1', 'cp_ps2', 'ki_ps2', 'MSW', '-append');
    %% %根值校验【Roots Check】
    if strcmp(szCheck, 'Yes') && ~strcmp(szChkCho, 'Never')
        [cp_3, ki_3] = RootsCheck(f, MPM, cp_0, ki_0);
        cp_ps3(n,1:1:length(cp_3)) = cp_3; ki_ps3(n, 1:1:length(ki_3)) = ki_3;
        save('data.mat', 'cp_ps3', 'ki_ps3', '-append');
    end
    %% %后处理【Post Process】
    disp('*******************************************************');
    diary off; n = n+1;
end
disp('done!');

% --- Executes during object creation, after setting all properties.
function Run_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function nMode_cs_Callback(hObject, eventdata, handles)
% hObject    handle to nMode_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nMode_cs as text
%        str2double(get(hObject,'String')) returns contents of nMode_cs as a double


% --- Executes during object creation, after setting all properties.
function nMode_cs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nMode_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in szMode.
function szMode_Callback(hObject, eventdata, handles)
% hObject    handle to szMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns szMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from szMode


% --- Executes during object creation, after setting all properties.
function szMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to szMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in szBC.
function szBC_Callback(hObject, eventdata, handles)
% hObject    handle to szBC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns szBC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from szBC


% --- Executes during object creation, after setting all properties.
function szBC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to szBC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in szFun.
function szFun_Callback(hObject, eventdata, handles)
% hObject    handle to szFun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns szFun contents as cell array
%        contents{get(hObject,'Value')} returns selected item from szFun


% --- Executes during object creation, after setting all properties.
function szFun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to szFun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c_outer_Callback(hObject, eventdata, handles)
% hObject    handle to c_outer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_outer as text
%        str2double(get(hObject,'String')) returns contents of c_outer as a double


% --- Executes during object creation, after setting all properties.
function c_outer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_outer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function b_middle_Callback(hObject, eventdata, handles)
% hObject    handle to b_middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of b_middle as text
%        str2double(get(hObject,'String')) returns contents of b_middle as a double


% --- Executes during object creation, after setting all properties.
function b_middle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_inner_Callback(hObject, eventdata, handles)
% hObject    handle to a_inner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_inner as text
%        str2double(get(hObject,'String')) returns contents of a_inner as a double


% --- Executes during object creation, after setting all properties.
function a_inner_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_inner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dvm_Callback(hObject, eventdata, handles)
% hObject    handle to dvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dvm as text
%        str2double(get(hObject,'String')) returns contents of dvm as a double


% --- Executes during object creation, after setting all properties.
function dvm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dem_Callback(hObject, eventdata, handles)
% hObject    handle to dem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dem as text
%        str2double(get(hObject,'String')) returns contents of dem as a double


% --- Executes during object creation, after setting all properties.
function dem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kur_Callback(hObject, eventdata, handles)
% hObject    handle to kur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kur as text
%        str2double(get(hObject,'String')) returns contents of kur as a double


% --- Executes during object creation, after setting all properties.
function kur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function err_Callback(hObject, eventdata, handles)
% hObject    handle to err (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of err as text
%        str2double(get(hObject,'String')) returns contents of err as a double


% --- Executes during object creation, after setting all properties.
function err_CreateFcn(hObject, eventdata, handles)
% hObject    handle to err (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kib_Callback(hObject, eventdata, handles)
% hObject    handle to kib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kib as text
%        str2double(get(hObject,'String')) returns contents of kib as a double


% --- Executes during object creation, after setting all properties.
function kib_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dki_Callback(hObject, eventdata, handles)
% hObject    handle to dki (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dki as text
%        str2double(get(hObject,'String')) returns contents of dki as a double


% --- Executes during object creation, after setting all properties.
function dki_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dki (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kia_Callback(hObject, eventdata, handles)
% hObject    handle to kia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kia as text
%        str2double(get(hObject,'String')) returns contents of kia as a double


% --- Executes during object creation, after setting all properties.
function kia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpb_Callback(hObject, eventdata, handles)
% hObject    handle to cpb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpb as text
%        str2double(get(hObject,'String')) returns contents of cpb as a double


% --- Executes during object creation, after setting all properties.
function cpb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dcp_Callback(hObject, eventdata, handles)
% hObject    handle to dcp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dcp as text
%        str2double(get(hObject,'String')) returns contents of dcp as a double


% --- Executes during object creation, after setting all properties.
function dcp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dcp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cpa_Callback(hObject, eventdata, handles)
% hObject    handle to cpa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cpa as text
%        str2double(get(hObject,'String')) returns contents of cpa as a double


% --- Executes during object creation, after setting all properties.
function cpa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cpa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fb_Callback(hObject, eventdata, handles)
% hObject    handle to fb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fb as text
%        str2double(get(hObject,'String')) returns contents of fb as a double


% --- Executes during object creation, after setting all properties.
function fb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function df_Callback(hObject, eventdata, handles)
% hObject    handle to df (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of df as text
%        str2double(get(hObject,'String')) returns contents of df as a double


% --- Executes during object creation, after setting all properties.
function df_CreateFcn(hObject, eventdata, handles)
% hObject    handle to df (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fa_Callback(hObject, eventdata, handles)
% hObject    handle to fa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fa as text
%        str2double(get(hObject,'String')) returns contents of fa as a double


% --- Executes during object creation, after setting all properties.
function fa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c2_Callback(hObject, eventdata, handles)
% hObject    handle to c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c2 as text
%        str2double(get(hObject,'String')) returns contents of c2 as a double


% --- Executes during object creation, after setting all properties.
function c2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function row2_Callback(hObject, eventdata, handles)
% hObject    handle to row2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of row2 as text
%        str2double(get(hObject,'String')) returns contents of row2 as a double


% --- Executes during object creation, after setting all properties.
function row2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to row2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c1_Callback(hObject, eventdata, handles)
% hObject    handle to c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1 as text
%        str2double(get(hObject,'String')) returns contents of c1 as a double


% --- Executes during object creation, after setting all properties.
function c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function row1_Callback(hObject, eventdata, handles)
% hObject    handle to row1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of row1 as text
%        str2double(get(hObject,'String')) returns contents of row1 as a double


% --- Executes during object creation, after setting all properties.
function row1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to row1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmavm_Callback(hObject, eventdata, handles)
% hObject    handle to sigmavm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmavm as text
%        str2double(get(hObject,'String')) returns contents of sigmavm as a double


% --- Executes during object creation, after setting all properties.
function sigmavm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmavm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ytavm_Callback(hObject, eventdata, handles)
% hObject    handle to ytavm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ytavm as text
%        str2double(get(hObject,'String')) returns contents of ytavm as a double


% --- Executes during object creation, after setting all properties.
function ytavm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ytavm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Evm0_Callback(hObject, eventdata, handles)
% hObject    handle to Evm0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Evm0 as text
%        str2double(get(hObject,'String')) returns contents of Evm0 as a double


% --- Executes during object creation, after setting all properties.
function Evm0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Evm0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rowvm_Callback(hObject, eventdata, handles)
% hObject    handle to rowvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rowvm as text
%        str2double(get(hObject,'String')) returns contents of rowvm as a double


% --- Executes during object creation, after setting all properties.
function rowvm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rowvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigmaem_Callback(hObject, eventdata, handles)
% hObject    handle to sigmaem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigmaem as text
%        str2double(get(hObject,'String')) returns contents of sigmaem as a double


% --- Executes during object creation, after setting all properties.
function sigmaem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ytaem_Callback(hObject, eventdata, handles)
% hObject    handle to ytaem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ytaem as text
%        str2double(get(hObject,'String')) returns contents of ytaem as a double


% --- Executes during object creation, after setting all properties.
function ytaem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ytaem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Eem0_Callback(hObject, eventdata, handles)
% hObject    handle to Eem0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Eem0 as text
%        str2double(get(hObject,'String')) returns contents of Eem0 as a double


% --- Executes during object creation, after setting all properties.
function Eem0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Eem0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rowem_Callback(hObject, eventdata, handles)
% hObject    handle to rowem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rowem as text
%        str2double(get(hObject,'String')) returns contents of rowem as a double


% --- Executes during object creation, after setting all properties.
function rowem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rowem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when uipanel8 is resized.
function uipanel8_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
