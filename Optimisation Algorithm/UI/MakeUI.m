function [SelectedToolbox, num_modes, max_ops, trunc, output_loss] = MakeUI (FullToolbox, ToolboxTypes, default_num_modes, default_max_operations, default_trunc, default_output_loss)

%=================== Rosanna Nichols =======================%
%========== https://arxiv.org/abs/1812.01032 ===============%

close all

%% SETUP

%set sizes and spacings
Checkbox_spacing = 30;
Setting_title_spacing = 25;
Setting_spacing = 22;
Panel_left = 20;
Panel_space_buffer = 30;
Panel_height_buffer = 40;
Panel_width = 270;
Panel_gap = 40;
Panel_horizontal_spacing = Panel_width + Panel_gap;
Panel_dist_from_bottom = 100;
Panel_dist_from_top = 50;
Button_margin = 10;
Button_height = 50;
Button_width = 180;
Button_bottom = (Panel_dist_from_bottom - Button_height)/2; 
Figure_width = 4*Panel_width + 3*Panel_gap + 2*Panel_left;
Other_settings_height = 7*Checkbox_spacing + 4*Setting_spacing;

%open figure window
ui_window = figure(1);

ui_window.Name = 'Select Toolbox';

% initialise cell arrays to hold ui elements we will need to access the
% responses of.
Checkboxes = cell(1,3);
MinMaxEdits = cell(1,3);
FixedSettingEdits = cell(1,3);

% find maximum height needed for all the elements of a panel
heights = zeros(4,1);
for type_index = 1:3
    no_modules = size(FullToolbox{type_index},1);
    heights(type_index) = heights(type_index) + no_modules * Checkbox_spacing;
    for module_index = 1:no_modules
        no_settings = FullToolbox{type_index}{module_index, 3};
        no_fixed_settings = FullToolbox{type_index}{module_index, 7};
        heights(type_index) = heights(type_index) + no_settings * (Setting_title_spacing + Setting_spacing) + no_fixed_settings * Setting_spacing;
    end
end
heights(4) = Other_settings_height;

max_height = max(heights) + Panel_space_buffer;

%Resize figure
Figure_height = max_height + Panel_dist_from_bottom + Panel_height_buffer + Panel_dist_from_top;
ui_window.Units = 'pixels';
ui_window.Position = [100, 140, Figure_width, Figure_height];


%Figure title (text)
Figure_title_width = 500;
Figure_title_height = 30;
Figure_title_left = 2*Panel_width + 1.5*Panel_gap + Panel_left - Figure_title_width/2;
Figure_title_bottom = Figure_height - (Panel_dist_from_top + Figure_title_height)/2;
uicontrol('Style','text',...
    'String', 'Toolbox Selection',...
    'HorizontalAlignment','center',...
    'Position', [Figure_title_left,Figure_title_bottom,Figure_title_width,Figure_title_height],...
    'FontSize', 20);


%% MODULE SELECTION PANELS

% Loop over module type, modules within those types and settings within
% those modules and create UI elements for each
for type_index = 1:3
    
    %Panel to contain modules of current module type
    Panel = uipanel('Title',ToolboxTypes{type_index},'FontSize',14,...
        'Units', 'pixels',...
        'Position',[Panel_left + Panel_horizontal_spacing*(type_index-1), Panel_dist_from_bottom, Panel_width, max_height + Panel_height_buffer]);
    
    %find no. of modules for this type, and initialise ui element arrays
    no_modules = size(FullToolbox{type_index},1);
    Checkboxes{type_index} = gobjects(no_modules,1);
    MinMaxEdits{type_index} = cell(no_modules,1);
    FixedSettingEdits{type_index} = cell(no_modules,1);
    
    %set initial distance from the bottom of the panel
    current_dist_from_bottom = max_height;
    
    %loop over modules within type
    for module_index = 1:no_modules
        
        %create checkbox
        current_dist_from_bottom = current_dist_from_bottom - Checkbox_spacing;
        Checkboxes{type_index}(module_index) = uicontrol(Panel,...
            'Style','checkbox',...
            'String', FullToolbox{type_index}{module_index, 1},...
            'Value', 1,...
            'Position', [10, current_dist_from_bottom, 240, 30],...
            'FontSize', 12);
        
        %find no. of settings for this module, and initialise ui element arrays
        no_settings = FullToolbox{type_index}{module_index, 3};
        MinMaxEdits{type_index}{module_index} = gobjects(no_settings, 2);
        
        %loop over settings of this module
        for setting_index = 1:no_settings
            
            %Setting title (text)
            current_dist_from_bottom = current_dist_from_bottom - Setting_title_spacing;
            uicontrol(Panel,'Style','text',...
                'String', FullToolbox{type_index}{module_index, 4}{setting_index},...
                'HorizontalAlignment','left',...
                'Position', [15,current_dist_from_bottom,100,20],...
                'FontSize', 10,...
                'FontWeight','bold');
            
            % "Min." (text)
            current_dist_from_bottom = current_dist_from_bottom - Setting_spacing;
            uicontrol(Panel,'Style','text',...
                'String', 'Min.',...
                'HorizontalAlignment','left',...
                'Position', [15,current_dist_from_bottom,35,15],...
                'FontSize', 10);
            
            % Editable field containing default minimum for that setting
            MinMaxEdits{type_index}{module_index}(setting_index, 1) = uicontrol(Panel, 'Style', 'edit',...
                'String', FullToolbox{type_index}{module_index, 5}(setting_index, 1),...
                'Position', [50 ,current_dist_from_bottom,41,22],...
                'FontSize', 10);
            
            % "Max." (text)
            uicontrol(Panel,'Style','text',...
                'String', 'Max.',...
                'HorizontalAlignment','left',...
                'Position', [105,current_dist_from_bottom,35,15],...
                'FontSize', 10);
            
            % Editable field containing default maximum for that setting
            MinMaxEdits{type_index}{module_index}(setting_index, 2) = uicontrol(Panel, 'Style', 'edit',...
                'String', FullToolbox{type_index}{module_index, 5}(setting_index, 2),...
                'Position', [145,current_dist_from_bottom,41,22],...
                'FontSize', 10);
            
        end
        
        %find no. of fixed settings for this module, and initialise ui element arrays
        no_fixed_settings = FullToolbox{type_index}{module_index, 7};
        FixedSettingEdits{type_index}{module_index} = gobjects(no_fixed_settings, 1);
        
        %loop over settings of this module
        for setting_index = 1:no_fixed_settings
        
            %Setting title (text)
            current_dist_from_bottom = current_dist_from_bottom - Setting_title_spacing;
            uicontrol(Panel,'Style','text',...
                'String', FullToolbox{type_index}{module_index, 8}{setting_index},...
                'HorizontalAlignment','left',...
                'Position', [15,current_dist_from_bottom,140,20],...
                'FontSize', 10,...
                'FontWeight','bold');
            
            % Editable field containing default value for that setting
            FixedSettingEdits{type_index}{module_index}(setting_index) = uicontrol(Panel, 'Style', 'edit',...
                'String', FullToolbox{type_index}{module_index, 9}(setting_index),...
                'Position', [145 ,current_dist_from_bottom,41,22],...
                'FontSize', 10);
            
        end
    end
end

%% OTHER SETTINGS PANEL

% Create number of modes, maximum number of operations, and truncation panel

Panel = uipanel('Title', 'Other options','FontSize',14,...
        'Units', 'pixels',...
        'Position',[Panel_left + 3*Panel_horizontal_spacing, Panel_dist_from_bottom, Panel_width, max_height + Panel_height_buffer]);

current_dist_from_bottom = max_height;
    
%Number of modes - title
current_dist_from_bottom = current_dist_from_bottom - Checkbox_spacing;
uicontrol(Panel,'Style','text',...
    'String', 'Number of modes',...
    'HorizontalAlignment','left',...
    'Position', [15,current_dist_from_bottom,200,20],...
    'FontSize', 12);

% Number of modes - editable
current_dist_from_bottom = current_dist_from_bottom - Setting_spacing;
NumModesEdit = uicontrol(Panel, 'Style', 'edit',...
    'String', default_num_modes,...
    'Position', [Panel_width/2 - 20.5 ,current_dist_from_bottom,41,22],...
    'FontSize', 10);

% Max. no. operations - title
current_dist_from_bottom = current_dist_from_bottom - 2*Checkbox_spacing;
uicontrol(Panel,'Style','text',...
    'String', 'Max. number of operations',...
    'HorizontalAlignment','left',...
    'Position', [15,current_dist_from_bottom,200,20],...
    'FontSize', 12);

% Max. no. operations - editable
current_dist_from_bottom = current_dist_from_bottom - Setting_spacing;
MaxOpsEdit = uicontrol(Panel, 'Style', 'edit',...
    'String', default_max_operations,...
    'Position', [Panel_width/2 - 20.5 ,current_dist_from_bottom,41,22],...
    'FontSize', 10);

% Truncation - title
current_dist_from_bottom = current_dist_from_bottom - 2*Checkbox_spacing;
uicontrol(Panel,'Style','text',...
    'String', 'Max. Truncation',...
    'HorizontalAlignment','left',...
    'Position', [15,current_dist_from_bottom,200,20],...
    'FontSize', 12);

current_dist_from_bottom = current_dist_from_bottom - Setting_spacing;
% Truncation - editable
TruncEdit = uicontrol(Panel, 'Style', 'edit',...
    'String', default_trunc,...
    'Position', [Panel_width/2 - 20.5, current_dist_from_bottom,41,22],...
    'FontSize', 10);

% Output loss - title
current_dist_from_bottom = current_dist_from_bottom - 2*Checkbox_spacing;
uicontrol(Panel,'Style','text',...
    'String', 'Loss on output state',...
    'HorizontalAlignment','left',...
    'Position', [15,current_dist_from_bottom,200,20],...
    'FontSize', 12);

current_dist_from_bottom = current_dist_from_bottom - Setting_spacing;
% Output loss - editable
OutputLossEdit = uicontrol(Panel, 'Style', 'edit',...
    'String', default_output_loss,...
    'Position', [Panel_width/2 - 20.5, current_dist_from_bottom,41,22],...
    'FontSize', 10);

%% BUTTONS

% Create "Load settings" button
uicontrol('Style','pushbutton','String','Load settings',...
           'Position',[Figure_width - (4*Button_margin + 3*Button_width), Button_bottom, Button_width, Button_height],...
           'Callback',  @LoadButton_Callback,...
           'FontSize', 20);

% Create "Save settings" button
uicontrol('Style','pushbutton','String','Save settings',...
           'Position',[Figure_width - (3*Button_margin + 2*Button_width), Button_bottom, Button_width, Button_height],...
           'Callback',  @SaveButton_Callback,...
           'FontSize', 20);

% Create "Next" button
uicontrol('Style','pushbutton','String','Next',...
           'Position',[Figure_width - (2*Button_margin + Button_width), Button_bottom, Button_width, Button_height],...
           'Callback',  @NextButton_Callback,...
           'FontSize', 20);


SelectedToolbox = cell(1,3);
% wait for figure to close (i.e. for go button to be clicked and input data
% to be collected) before ending the function
uiwait

%% CALLBACK FUNCTIONS

    % when "Save settings" button is clicked, put user input in selected toolbox
    function SaveButton_Callback (source, eventdata)
        [SelectedToolbox, num_modes, max_ops, trunc, output_loss] = GetInput(Checkboxes, MinMaxEdits, FixedSettingEdits, FullToolbox, NumModesEdit, MaxOpsEdit, TruncEdit, OutputLossEdit);
        uisave({'SelectedToolbox', 'num_modes', 'max_ops', 'trunc', 'output_loss'},'saved_toolbox_settings.mat')
    end

    % when "Load settings" button is clicked, put user input in selected toolbox
    function LoadButton_Callback (source, eventdata)
        [file, path] = uigetfile('*.mat', 'Select settings file');
        if file ~= 0
            loaded_vars = load(fullfile(path,file));
            SelectedToolbox = loaded_vars.SelectedToolbox;
            num_modes = loaded_vars.num_modes;
            max_ops = loaded_vars.max_ops;
            trunc = loaded_vars.trunc;
            output_loss = loaded_vars.output_loss;
            close(gcbf)
        end
    end

    % when "Next" button is clicked, put user input in selected toolbox and
    % close figure
    function NextButton_Callback (source, eventdata)
        [SelectedToolbox, num_modes, max_ops, trunc, output_loss] = GetInput(Checkboxes, MinMaxEdits, FixedSettingEdits, FullToolbox, NumModesEdit, MaxOpsEdit, TruncEdit, OutputLossEdit);
        close(gcbf)
    end


end