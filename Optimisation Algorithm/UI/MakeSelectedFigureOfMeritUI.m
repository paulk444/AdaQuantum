function MakeSelectedFigureOfMeritUI (FiguresOfMerit, use_heralding_probability)

%=================== Rosanna Nichols =======================%
%========== https://arxiv.org/abs/1812.01032 ===============%

    %set sizes and spacings
    Checkbox_spacing = 30;
    Setting_spacing = 22;
    Panel_left = 20;
    Panel_space_buffer = 30;
    Panel_height_buffer = 40;
    Panel_width = 270;
    Panel_gap = 40;
    Panel_dist_from_bottom = 100;
    Panel_dist_from_top = 50;
    Figure_width = 3*Panel_width + 2*Panel_gap + 2*Panel_left;
    no_aug_functions = size(FiguresOfMerit{2},1);
    Aug_panel_height = 320 + no_aug_functions * (Checkbox_spacing + 2 * Setting_spacing) + 27;
    Setting_horizontal_spacing = 120;

    %open figure window
    ui_window = figure(2);
    ui_window.Name = 'Chosen Figure of merit';

    %Resize figure
    Figure_height = Aug_panel_height + Panel_dist_from_bottom + Panel_height_buffer + Panel_dist_from_top;
    ui_window.Units = 'pixels';
    ui_window.Position = [100, 140, Figure_width, Figure_height];

    %Figure title (text)
    Figure_title_width = 500;
    Figure_title_height = 30;
    Figure_title_left = 1.5*Panel_width + Panel_gap + Panel_left - Figure_title_width/2;
    Figure_title_bottom = Figure_height - (Panel_dist_from_top + Figure_title_height)/2;
    uicontrol('Style','text',...
        'String', 'Figure of Merit Selection',...
        'HorizontalAlignment','center',...
        'Position', [Figure_title_left,Figure_title_bottom,Figure_title_width,Figure_title_height],...
        'FontSize', 20);

    %%Radio buttons for figure of merit selection
    no_FoMs = size(FiguresOfMerit{1},1);
    max_height_FoM = no_FoMs * Checkbox_spacing + Panel_height_buffer;
    ButtonGroup_height = max_height_FoM + Panel_height_buffer;
    current_dist_from_bottom = Figure_height - Panel_dist_from_top - ButtonGroup_height - Panel_height_buffer;
    ButtonGroup = uibuttongroup('Title', 'Select Figure of Merit',...
        'FontSize', 14,...
        'Units', 'pixels',...
        'Position', [Panel_left, current_dist_from_bottom, Panel_width, ButtonGroup_height]);
    FigureOfMeritButtons = gobjects(no_FoMs, 1);
    inner_current_dist_from_bottom = max_height_FoM;
    for FoM_index = 1:no_FoMs

        inner_current_dist_from_bottom = inner_current_dist_from_bottom - Checkbox_spacing;

        FigureOfMeritButtons(FoM_index) = uicontrol(ButtonGroup,...
                'Style','radiobutton',...
                'String', FiguresOfMerit{1}{FoM_index, 1},...
                'Position', [10, inner_current_dist_from_bottom, 240, 30],...
                'FontSize', 12,...
                'Enable', 'off');

    end

    %% Heralding prob

    Heralding_prob_panel_height = 2*Checkbox_spacing + Panel_height_buffer;
    current_dist_from_bottom = current_dist_from_bottom - Heralding_prob_panel_height - Panel_height_buffer;
    %Panel to contain modules of current module type
    Heralding_prob_panel = uipanel('Title','Heralding probability','FontSize',14,...
        'Units', 'pixels',...
        'Position',[Panel_left, current_dist_from_bottom, Panel_width, Heralding_prob_panel_height]);

    uicontrol(Heralding_prob_panel,...
                'Style','checkbox',...
                'String','Use heralding probability',...
                'Value', use_heralding_probability,...
                'Position', [10, Checkbox_spacing, 240, Checkbox_spacing],...
                'FontSize', 12,...
                'Enable', 'off');

    %% Augmenting function

    %Panel to contain augmenting functions
    Aug_panel = uibuttongroup('Title','Augmenting function','FontSize',14,...
        'Units', 'pixels',...
        'Position',[Panel_left + Panel_width + Panel_gap, Panel_dist_from_bottom, 2*Panel_width, Aug_panel_height + Panel_space_buffer]);

    AugmentingFunctionButtons = gobjects(no_aug_functions, 1);
    inner_current_dist_from_bottom = Aug_panel_height;

    AugSettingEdits = cell(no_aug_functions,1);

    for Aug_index = 1:no_aug_functions

        inner_current_dist_from_bottom = inner_current_dist_from_bottom - Checkbox_spacing;

        AugmentingFunctionButtons(Aug_index) = uicontrol(Aug_panel,...
            'Style','radiobutton',...
            'String', FiguresOfMerit{2}{Aug_index, 1},...
            'Position', [10, inner_current_dist_from_bottom, 240, 30],...
            'FontSize', 12,...
            'Enable', 'off');

        num_settings = FiguresOfMerit{2}{Aug_index, 2};
        inner_current_dist_from_bottom = inner_current_dist_from_bottom - Setting_spacing + 5;

        AugSettingEdits{Aug_index} = gobjects(num_settings, 1);

        for setting_index = 1:num_settings

            uicontrol(Aug_panel,'Style','text',...
                'String', FiguresOfMerit{2}{Aug_index, 3}{setting_index},...
                'HorizontalAlignment','left',...
                'Position', [15 + (setting_index - 1) * Setting_horizontal_spacing, inner_current_dist_from_bottom, Setting_horizontal_spacing, 15],...
                'FontSize', 10);
        end

        inner_current_dist_from_bottom = inner_current_dist_from_bottom - Setting_spacing - 5;

        for setting_index = 1:num_settings

            AugSettingEdits{Aug_index}(setting_index) = uicontrol(Aug_panel,'Style','edit',...
                'String', FiguresOfMerit{2}{Aug_index, 4}(setting_index),...
                'HorizontalAlignment','left',...
                'Position', [15 + (setting_index - 1) * Setting_horizontal_spacing, inner_current_dist_from_bottom, 41, 22],...
                'FontSize', 10,...
                'Enable', 'off');
        end

    end


    % augmenting function preview axes
    axes('Parent', Aug_panel, 'Units', 'pixels', 'Position', [55,55,350,250], 'XLim', [0,5], 'YLim', [0,1.05]);
    xlabel('nbar')
    ylabel('Augmenting function value')

    PreviewCallback();

    function PreviewCallback (~,~)
        
        for i = 1:no_aug_functions
            if AugmentingFunctionButtons(i).Value

                PreviewFunction = FiguresOfMerit{2}{i, 5};

            end
        end
        
        nbar_vals = 0:0.01:5;
        PreviewFunctionVals = zeros(length(nbar_vals));
        for k = 1:length(nbar_vals)
            PreviewFunctionVals(k) = PreviewFunction(nbar_vals(k));
        end
        axes('Parent', Aug_panel, 'Units', 'pixels', 'Position', [55,55,350,250], 'XLim', [0,5], 'YLim', [0,1.05]);
        xlabel('nbar')
        ylabel('Augmenting function value')
        % hold so the plot uses these axes and doesn't draw new ones
        hold on
        plot(nbar_vals, PreviewFunctionVals)
        
    end


end