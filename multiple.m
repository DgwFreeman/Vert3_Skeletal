function varargout = multiple(varargin)
% MULTIPLE MATLAB code for multiple.fig
%      MULTIPLE, by itself, creates a new MULTIPLE or raises the existing
%      singleton*.
%
%      H = MULTIPLE returns the handle to a new MULTIPLE or the handle to
%      the existing singleton*.
%
%      MULTIPLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTIPLE.M with the given input arguments.
%
%      MULTIPLE('Property','Value',...) creates a new MULTIPLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multiple_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multiple_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help multiple

% Last Modified by GUIDE v2.5 01-May-2014 13:43:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multiple_OpeningFcn, ...
                   'gui_OutputFcn',  @multiple_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before multiple is made visible.
function multiple_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multiple (see VARARGIN)

% Choose default command line output for multiple
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes multiple wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = multiple_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function AppendToListbox(handle, string)
% Put 'string' on the end of 'handle's values
contents = cellstr(get(handle, 'String'));
new_entry = {string};
set(handle, 'String', char([contents; new_entry]));

function vec = ListboxToVector(listbox)
contents = cellstr(get(listbox, 'String'));
indices = get(listbox, 'Value');
% had some trouble finding something that works if indices = []
% but this is simple enough
vec = cellfun(@str2double, contents(indices));
% need to return a row vector for matlab's for, i think
vec = vec';

function name = OutdirName(~, ~, ~, ~, ~, ~)
% This function is the worst
% Query the directory to find the next available set # and use that
% FIXME replace with name thaat uses parameters, or better, a non-directory
% organization
existent = dir(['DataFiles' filesep 'Multiple*']);
if (isempty(existent))
    lastnum = 0;
else
    % used to use existent(end) here but 10 < 9 in alphabetic sort
    lastnum = max(arrayfun(@(v) str2double(v.name(length('Multiple')+1:end)), existent));
end
name = ['Multiple' num2str(lastnum+1)];
    
% --- Executes on selection change in TnDensityRandom.
function TnDensityRandom_Callback(~, ~, ~)
% hObject    handle to TnDensityRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TnDensityRandom contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TnDensityRandom


% --- Executes during object creation, after setting all properties.
function TnDensityRandom_CreateFcn(hObject, ~, ~)
% hObject    handle to TnDensityRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TnDensityUniform.
function TnDensityUniform_Callback(~, ~, ~)
% hObject    handle to TnDensityUniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TnDensityUniform contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TnDensityUniform


% --- Executes during object creation, after setting all properties.
function TnDensityUniform_CreateFcn(hObject, ~, ~)
% hObject    handle to TnDensityUniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TnDensityRandomEdit_Callback(~, ~, ~)
% hObject    handle to TnDensityRandomEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TnDensityRandomEdit as text
%        str2double(get(hObject,'String')) returns contents of TnDensityRandomEdit as a double


% --- Executes during object creation, after setting all properties.
function TnDensityRandomEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to TnDensityRandomEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TnDensityRandomAdd.
function TnDensityRandomAdd_Callback(hObject, eventdata, handles)
% hObject    handle to TnDensityRandomAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
adding = get(handles.TnDensityRandomEdit, 'String');
addn = str2double(adding);
if (isnan(addn) || addn < 0 || addn > 1) % NaN should compare false but let's be explicit
    warndlg(sprintf('Densities should be numbers between 0 and 1, not %s.', adding), 'Bad input');
else
    AppendToListbox(handles.TnDensityRandom, adding);
end

function TnDensityUniformEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TnDensityUniformEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TnDensityUniformEdit as text
%        str2double(get(hObject,'String')) returns contents of TnDensityUniformEdit as a double


% --- Executes during object creation, after setting all properties.
function TnDensityUniformEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TnDensityUniformEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TnDensityUniformAdd.
function TnDensityUniformAdd_Callback(hObject, eventdata, handles)
% hObject    handle to TnDensityUniformAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
adding = get(handles.TnDensityUniformEdit, 'String');
addn = str2double(adding);
if (isnan(addn) || addn < 0 || addn > 1)
    warndlg(sprintf('Densities should be numbers between 0 and 1, not %s.', adding), 'Bad input');
else
    AppendToListbox(handles.TnDensityUniform, adding);
end

% --- Executes on selection change in XBDensityRandom.
function XBDensityRandom_Callback(hObject, eventdata, handles)
% hObject    handle to XBDensityRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns XBDensityRandom contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XBDensityRandom


% --- Executes during object creation, after setting all properties.
function XBDensityRandom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XBDensityRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in XBDensityUniform.
function XBDensityUniform_Callback(hObject, eventdata, handles)
% hObject    handle to XBDensityUniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns XBDensityUniform contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XBDensityUniform


% --- Executes during object creation, after setting all properties.
function XBDensityUniform_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XBDensityUniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XBDensityRandomEdit_Callback(hObject, eventdata, handles)
% hObject    handle to XBDensityRandomEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XBDensityRandomEdit as text
%        str2double(get(hObject,'String')) returns contents of XBDensityRandomEdit as a double


% --- Executes during object creation, after setting all properties.
function XBDensityRandomEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XBDensityRandomEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in XBDensityRandomAdd.
function XBDensityRandomAdd_Callback(~, ~, handles)
% hObject    handle to XBDensityRandomAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
adding = get(handles.XBDensityRandomEdit, 'String');
addn = str2double(adding);
if (isnan(addn) || addn < 0 || addn > 1)
    warndlg(sprintf('Densities should be numbers between 0 and 1, not %s.', adding), 'Bad input');
else
    AppendToListbox(handles.XBDensityRandom, adding);
end

function XBDensityUniformEdit_Callback(~, ~, ~)
% hObject    handle to XBDensityUniformEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XBDensityUniformEdit as text
%        str2double(get(hObject,'String')) returns contents of XBDensityUniformEdit as a double


% --- Executes during object creation, after setting all properties.
function XBDensityUniformEdit_CreateFcn(hObject, ~, ~)
% hObject    handle to XBDensityUniformEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in XBDensityUniformAdd.
function XBDensityUniformAdd_Callback(hObject, eventdata, handles)
% hObject    handle to XBDensityUniformAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
adding = get(handles.XBDensityUniformEdit, 'String');
addn = str2double(adding);
if (isnan(addn) || addn < 0 || addn > 1)
    warndlg(sprintf('Densities should be numbers between 0 and 1, not %s.', adding), 'Bad input');
else
    AppendToListbox(handles.XBDensityUniform, adding);
end


% --- Executes on selection change in pCa.
function pCa_Callback(~, ~, ~)
% hObject    handle to pCa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pCa contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pCa


% --- Executes during object creation, after setting all properties.
function pCa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pCa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in kxscaler.
function kxscaler_Callback(hObject, eventdata, handles)
% hObject    handle to kxscaler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns kxscaler contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kxscaler


% --- Executes during object creation, after setting all properties.
function kxscaler_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kxscaler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HalfSL_Callback(hObject, eventdata, handles)
% hObject    handle to HalfSL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HalfSL as text
%        str2double(get(hObject,'String')) returns contents of HalfSL as a double


% --- Executes during object creation, after setting all properties.
function HalfSL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HalfSL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LTitin_Callback(~, ~, ~)
% hObject    handle to LTitin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LTitin as text
%        str2double(get(hObject,'String')) returns contents of LTitin as a double


% --- Executes during object creation, after setting all properties.
function LTitin_CreateFcn(hObject, ~, ~)
% hObject    handle to LTitin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GoButton.
function GoButton_Callback(hObject, ~, handles)
% hObject    handle to GoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

TnDensitiesUniform = ListboxToVector(handles.TnDensityUniform);
TnDensitiesRandom = ListboxToVector(handles.TnDensityRandom);
XBDensitiesUniform = ListboxToVector(handles.XBDensityUniform);
XBDensitiesRandom = ListboxToVector(handles.XBDensityRandom);
pCaV = ListboxToVector(handles.pCa);
kxscaler = ListboxToVector(handles.kxscaler);
HalfSL = str2double(get(handles.HalfSL, 'String'));
LTitin = str2double(get(handles.LTitin, 'String'));

NumRuns = str2double(get(handles.NumRuns, 'String'));

if(isempty(TnDensitiesUniform) && isempty(TnDensitiesRandom))
    warndlg('Please select at least one troponin density.', 'Insufficient input');
elseif (isempty(XBDensitiesUniform) && isempty(XBDensitiesRandom))
    warndlg('Please select at least one cross-bridge density.', 'Insufficient input');
elseif (isempty(pCaV))
    warndlg('Please select at least one calcium level.', 'Insufficient input');
elseif (isempty(kxscaler))
    warndlg('Please select at least one cross-bridge stiffness scaler.', 'Insufficient input');
else
    %% FIXME.
    init_params;
    for TnKOType=[0,1]
        filaments.TnKOType = TnKOType;
        if TnKOType == 0
            TnDensity = TnDensitiesRandom;
        else
            TnDensity = TnDensitiesUniform;
        end
        for TnKO = TnDensity
            filaments.TnFraction = TnKO;
            for XBKOType=[0,1]
                filaments.XBKOType = XBKOType;
                if XBKOType == 0
                    XBDensity = XBDensitiesRandom;
                else
                    XBDensity = XBDensitiesUniform;
                end
                for XBKO = XBDensity
                    filaments.XBFraction = XBKO;
                    for kxscaler = kxscaler
                        StiffScale.kxscaler = kxscaler;
                        OutDir = ['DataFiles' filesep OutdirName(TnKOType, TnDensity, XBKOType, XBDensity, kxscaler, HalfSL) filesep];
                        mkdir(OutDir);
                        for pCa_in = pCaV
                            % pCa is a parameter by itself, i.e. no struct
                            [Steps, Means, Vars, IndexThalf, Binder] = RunSeveral(NumRuns, DataParams, StartLength, pCa_in, StiffScale, filaments, knockout, coop, TFRateScale, tcparam);
                            WriteText(OutDir, pCa_in, DataParams.dt, Binder, Steps, Means, Vars, IndexThalf);
                        end
                    end
                end
            end
        end
    end
end

function NumRuns_Callback(hObject, eventdata, handles)
% hObject    handle to NumRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumRuns as text
%        str2double(get(hObject,'String')) returns contents of NumRuns as a double


% --- Executes during object creation, after setting all properties.
function NumRuns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumRuns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
