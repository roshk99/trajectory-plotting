function varargout = canal_visualization_gui(varargin)
% CANAL_VISUALIZATION_GUI MATLAB code for canal_visualization_gui.fig
%      CANAL_VISUALIZATION_GUI, by itself, creates a new CANAL_VISUALIZATION_GUI or raises the existing
%      singleton*.
%
%      H = CANAL_VISUALIZATION_GUI returns the handle to a new CANAL_VISUALIZATION_GUI or the handle to
%      the existing singleton*.
%
%      CANAL_VISUALIZATION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CANAL_VISUALIZATION_GUI.M with the given input arguments.
%
%      CANAL_VISUALIZATION_GUI('Property','Value',...) creates a new CANAL_VISUALIZATION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before canal_visualization_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to canal_visualization_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help canal_visualization_gui

% Last Modified by GUIDE v2.5 20-Jun-2016 14:06:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @canal_visualization_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @canal_visualization_gui_OutputFcn, ...
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


% --- Executes just before canal_visualization_gui is made visible.
function canal_visualization_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to canal_visualization_gui (see VARARGIN)

% Choose default command line output for canal_visualization_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes canal_visualization_gui wait for user response (see UIRESUME)
% uiwait(handles.canal_visualization_gui);


% --- Outputs from this function are returned to the command line.
function varargout = canal_visualization_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in dataset_dropdown.
function dataset_dropdown_Callback(hObject, eventdata, handles)
% hObject    handle to dataset_dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dataset_dropdown contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataset_dropdown


% --- Executes during object creation, after setting all properties.
function dataset_dropdown_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataset_dropdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);

%Get Dataset Number
set_to_run = get(handles.dataset_dropdown,'Value') - 1;

%Get Fit Type
if get(handles.circles_radio, 'Value')
    fit_type = 'circles';
elseif get(handles.ellipses_radio, 'Value')
    fit_type = 'ellipses';
end

%Get Trajectories
how_many = [];
if get(handles.trajectory1, 'Value')
    how_many = [how_many 1];
end
if get(handles.trajectory2, 'Value')
    how_many = [how_many 2];
end
if get(handles.trajectory3, 'Value')
    how_many = [how_many 3];
end
if get(handles.trajectory4, 'Value') && set_to_run < 4
    how_many = [how_many 4];
end
if get(handles.trajectory5, 'Value') && set_to_run == 3
    how_many = [how_many 5];
end
if isempty(how_many)
    switch set_to_run
        case 0
            how_many = 1:4;
        case 1
            how_many = 1:4;
        case 2
            how_many = 1:4;
        case 3
            how_many = 1:5;
        case 4
            how_many = 1:5;
    end
end

%Get Plot Method
if get(handles.circles_plot_radio, 'Value')
    plot_method = 'circles';
else
    plot_method = 'surface';
end

runFunction(set_to_run, fit_type, plot_method, how_many, ...
    handles);

% --- Executes on button press in trajectory1.
function trajectory1_Callback(hObject, eventdata, handles)
% hObject    handle to trajectory1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trajectory1


% --- Executes on button press in trajectory2.
function trajectory2_Callback(hObject, eventdata, handles)
% hObject    handle to trajectory2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trajectory2


% --- Executes on button press in trajectory3.
function trajectory3_Callback(hObject, eventdata, handles)
% hObject    handle to trajectory3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trajectory3


% --- Executes on button press in trajectory4.
function trajectory4_Callback(hObject, eventdata, handles)
% hObject    handle to trajectory4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trajectory4


% --- Executes on button press in trajectory5.
function trajectory5_Callback(hObject, eventdata, handles)
% hObject    handle to trajectory5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trajectory5
