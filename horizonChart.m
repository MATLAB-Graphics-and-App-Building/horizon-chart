classdef horizonChart < matlab.graphics.chartcontainer.ChartContainer  & ...
        matlab.graphics.chartcontainer.mixin.Legend

    % horizonChart creates a Horizon Chart x-axis data and
    % y-axis data for each slice
    % 
    % horizonChart(x, y) creates a horizon chart with x an y data and
    % NumBands = 2 i.e. in each of the slices, the data is divided into
    % 2 regions - one above the baseline and one below the baseline
    % 
    % horizonChart(x, y, numBands) creates a horizon chart with x and y data
    % with number of bands = numBands. Number of Bands is the number of
    % sections that the data is divided into.
    % 
    % horizonChart(__ , Name, Value) specifies additional options for the horizon chart
    % using one or more name-value pair arguments. Specify the options after all other
    % input arguments.
        
    properties
        XData (:, :) double = [];
        YData (:, :) double = [];
        Labels (1, :) string = [];
        NumBands double {mustBeGreaterThanOrEqual(NumBands, 1)} = 2;

        XLabel string = "";
        YLabel string = "";
        Title string = "";
    end

    properties(Dependent)
        Baseline double = 0;
        ColorAboveBaseline {validatecolor} = [0, 0, 1];
        ColorBelowBaseline {validatecolor} = [1, 0, 0];
    end
    
    properties(Access = private,Transient,NonCopyable)
        PatchObjectListPerSlice (:, :)  = [];
        AxisHandlePerSlice (1, :) = [];
        SegRanges (1, :) double = [];
        NumSlices;
        ColorMap (:, 3) double; 
        
        
    end

    properties(Access = private)
        Baseline_I = NaN ;
        BaselineSetMethod string = "auto";
        
        ColorAboveBaseline_I = NaN;
        ColorAboveBaselineSetMethod = "auto";

        ColorBelowBaseline_I = NaN;
        ColorBelowBaselineSetMethod = "auto";
        
        ContainsMappingToolbox matlab.lang.OnOffSwitchState = 'on';
    end


    methods 
        function obj = horizonChart(varargin)
           % Intialize list of arguments
           args = varargin;
           leadingArgs = cell(0);

            if numel(args) >= 2 && isnumeric(args{1}) ...
                && isnumeric(args{2}) 
                   
                x = args{1};   
                y = args{2};

                if mod(numel(args), 2) == 1 && isnumeric(args{3})
                   % horizonChart(x, y, numBands)
                   numBands = args{3};
                   leadingArgs = [leadingArgs {'XData', x, 'YData', y , 'NumBands', numBands}];
                   args = args(4:end);
                else
                   % horizonChart(x, y)
                   leadingArgs = [leadingArgs {'XData', x, 'YData', y }];
                   args = args(3:end);
                end
            else
                if numel(args) < 2
                    error('Invalid Input Arguments. Too few inputs were provided. ')
                else
                    error('The first two arguments are not numeric and do not conform to XData and YData definition')
                end
            end
            
            % Combine positional arguments with name/value pairs
            args = [leadingArgs args];

            % Call superclass constructor method
            obj@matlab.graphics.chartcontainer.ChartContainer(args{:});

        end
    end
    methods(Access=protected)
        function setup(obj)
            if ~any(strcmp('Mapping Toolbox', {ver().Name})) 
                obj.ContainsMappingToolbox = 'off';
                warning("Mapping Toolbox is not installed. " + ...
                    "This may lead to degraded performance of the horizon chart. " + ...
                    "Install Mapping Toolbox for better performance")
            end
        end

        function update(obj)
            
            [obj.XData, obj.YData] = transformInputData(obj.XData, obj.YData) ;
            % Validate Inputs 
            validateInputs(obj.XData, obj.YData, obj.Labels);
            
            % If Grid Layout is already defined then clear the layout
            % during the update step
            children = obj.getLayout().Children;
            set(children, "Parent" , []);
            
            % Clear all patch objects 
            obj.PatchObjectListPerSlice = [];
            obj.AxisHandlePerSlice = [];            

            % Set GridLayout to be vertical layout
            obj.NumSlices = size(obj.YData, 2);
            obj.getLayout().GridSize = [obj.NumSlices 1];
            title(obj.getLayout(), obj.Title);
            xlabel(obj.getLayout(), obj.XLabel);
            ylabel(obj.getLayout(), obj.YLabel);
            
            % If user doesn't specify baseline we setup the baseline as the
            % median of the data. If the user specifies baseline, we adjust
            % the segment lengths to match the new baseline.
            if obj.BaselineSetMethod == "auto"
                obj.calculateSegmentsWithoutUserSetBaseline();
            else
                obj.calculateSegmentsWithUserSetBaseline();
            end 
            
            for slice = 1:obj.NumSlices
                % We obtain XData/ YData for each slice and calculate 
                % which band each data point belongs to  
                sliceXData = obj.XData(:, min(slice,size(obj.XData,2)))';
                sliceYData = obj.YData(:, slice)';

                binsPoint = binData(sliceYData, obj.SegRanges);

                PatchObjectList = [];
                order = [];
               
                % Get axis for the current tile and set all the properties
                ax = nexttile(obj.getLayout());
               
                % Disable data tips
                disableDefaultInteractivity(ax)

                % Specify labels only if the user has specified labels for
                % each slice
                if slice <= numel(obj.Labels) 
                    title(ax, obj.Labels(slice));
                end 

                y_min_slice = max(obj.YData(:));
                
                ax.XLim = [min(sliceXData), max(sliceXData)];
                ax.YTickLabel = [];
                ax.YTick = [];

                polygonOrder = gobjects(0);
                hold(ax, 'all')
                
                % Calculate color map for all the bands
                bgColor = get(ax, 'Color');
                order = ax.ColorOrder;
                
                if obj.ColorAboveBaselineSetMethod == "auto"
                    obj.ColorAboveBaseline_I = order(1, :); 
                end
                
                if obj.ColorBelowBaselineSetMethod == "auto"
                    obj.ColorBelowBaseline_I = order(2, :);
                end
               
                
                obj.CalculateColorMap(bgColor);


                % For each band we create a polygon that makes up the area
                % of the band. 
                for band = 1:obj.NumBands  
                    lower = obj.SegRanges(band);
                    upper = obj.SegRanges(band + 1);
                    color = obj.ColorMap(band, :);

                    % Calculate the vertices of the polygon depeding on
                    % whether it lies above the baseline/ or below the
                    % baseline 
                    [x_vertices, y_vertices] = generatePolygonPoints(sliceXData, sliceYData, lower, upper, lower >= obj.Baseline_I, obj.ContainsMappingToolbox);
                    
                    % Transform polygon and reflect it over the baseline
                    y_vertices = transformPolygon(y_vertices, lower, upper, obj.Baseline_I);
                                        
                    % Create the PatchObject for the band
                    PatchObject = patch(ax, x_vertices, y_vertices, color, 'DisplayName', num2str(lower) + " - " + num2str(upper));

                    
                    % If x_vertices/ y_vertices are empty then we need to
                    % create an empty patch object
                    if numel(x_vertices) == 0
                         PatchObject = patch(ax, NaN, NaN, color, 'DisplayName', num2str(lower) + " - " + num2str(upper));
                    else
                         % Find minimum all transformed y data in a particular slice
                         y_min_slice = min(y_min_slice, min(y_vertices(:)));
                    end  

                    % The bands that lie furthest from the baseline are
                    % displayed in the front. While the bands that are the
                    % closest to the baseline are displayed in the back
                    if lower >= obj.Baseline_I
                        polygonOrder = [PatchObject, polygonOrder];
                    else
                        polygonOrder = [polygonOrder, PatchObject];
                    end
                    
                    PatchObjectList = [PatchObjectList, PatchObject];
                end
                
                if y_min_slice ~= obj.Baseline_I
                    ax.YLim(1) = y_min_slice;
                end

                ax.Children = polygonOrder;
                hold(ax, 'off')

                obj.PatchObjectListPerSlice = [obj.PatchObjectListPerSlice; PatchObjectList];
                obj.AxisHandlePerSlice = [obj.AxisHandlePerSlice, ax];
            end  

            cbh = obj.buildColorBar(obj.AxisHandlePerSlice(end));
            cbh.Layout.Tile = 'east';
        end

        function CalculateColorMap(obj, backgroundColor)

            % The color of a band is decided by whether it lies above or
            % below the baseline. In case of bands that lie below the
            % baseline, the lower bands have a darker shade of
            % obj.colorsBelowBaseline. In case of bands that lie above the
            % baseline, the upper bands have a darker shade of
            % obj.colorsAboveBaseline
            nBandsBelowBaseline = sum(obj.SegRanges(2:end)<=obj.Baseline_I);
            nBandsAboveBaseline = obj.NumBands - nBandsBelowBaseline;
            
            % Calculate color gradient for the bands below the baseline
            alphas = fliplr(linspace(0.5, 1, nBandsBelowBaseline))';
            colorsBelowBaseline = alphas .* obj.ColorBelowBaseline_I + (1 - alphas) .* backgroundColor;

            % Calculate color gradient for the bands above the baseline
            alphas = linspace(0.5, 1, nBandsAboveBaseline)';
            colorsAboveBaseline = alphas .* obj.ColorAboveBaseline_I + (1 - alphas) .* backgroundColor;

            obj.ColorMap = [colorsBelowBaseline; colorsAboveBaseline];

        end
        
    end

    methods(Access = private)
        function calculateSegmentsWithoutUserSetBaseline(obj)
            % We divide the data into segments which contain equal amount of
            % data. For eg: If NumBands = 5, the first segment contains 20%
            % of the data. The second segment represents 20% - 40% of data
            % and so on
            obj.SegRanges = quantile(obj.YData(:), linspace(0, 1, obj.NumBands + 1));
            if mod(obj.NumBands, 2) == 0
                obj.Baseline_I = obj.SegRanges(obj.NumBands / 2 + 1);
            else
                obj.Baseline_I = obj.SegRanges((obj.NumBands + 1) / 2);
            end
        end

        function calculateSegmentsWithUserSetBaseline(obj)
            % In the first step, we divide data into segments using the
            % method proposed above.
            % Then we calculate segments below the baseline and above the
            % baseline. We accordingly divide the data below/above the
            % baseline using the newly found segments. 

            all_data = obj.YData(:);
            max_data = max(all_data);
            min_data = min(all_data);

            if obj.Baseline_I >= max_data
                nSegmentsAboveBaseline = 0;
                nSegmentsBelowBaseline = obj.NumBands;
            elseif obj.Baseline_I <= min_data
                nSegmentsAboveBaseline = obj.NumBands;
                nSegmentsBelowBaseline = 0;
            else
                segRanges = quantile(all_data, linspace(0, 1, obj.NumBands + 1));
                nSegmentsBelowBaseline = find(obj.Baseline_I >= segRanges, 1, 'last');

                if nSegmentsBelowBaseline == 0
                    nSegmentsBelowBaseline = nSegmentsBelowBaseline + 1;
                elseif nSegmentsBelowBaseline == obj.NumBands
                    nSegmentsBelowBaseline = nSegmentsBelowBaseline - 1;
                end

                nSegmentsAboveBaseline = obj.NumBands- nSegmentsBelowBaseline;
            end

            dataBelowBaseline = all_data(all_data < obj.Baseline_I);
            dataAboveBaseline = all_data(all_data > obj.Baseline_I);

            segRangesBelowBaseline = [];
            segRangesAboveBaseline = [];

            if nSegmentsBelowBaseline ~= 0 
                segRangesBelowBaseline = quantile(dataBelowBaseline, linspace(0, 1, nSegmentsBelowBaseline + 1));
            end 

            if nSegmentsAboveBaseline ~= 0
                segRangesAboveBaseline = quantile(dataAboveBaseline, linspace(0, 1, nSegmentsAboveBaseline + 1));
            end 

            obj.SegRanges = [segRangesBelowBaseline(1:nSegmentsBelowBaseline), obj.Baseline_I, segRangesAboveBaseline(2: nSegmentsAboveBaseline + 1)];
        end

        function cbh = buildColorBar(obj, ax)
            % Build a colorbar where the length of each color is proportional 
            % to the ratio of the length of each segment
            segLengths = obj.SegRanges(2:end) - obj.SegRanges(1: end - 1);
            lengthRatios = segLengths / sum(segLengths);
            lengthRatios = round(lengthRatios * 100);
            tLength = sum(lengthRatios);
            modifiedColorMap = colormap(repelem(obj.ColorMap, lengthRatios, 1));
            cbh = colorbar(ax);
            cumulLengthRatios = cumsum(lengthRatios) / tLength;
            cbh.Ticks = [0, cumulLengthRatios];
            cbh.TickLabels = num2cell(obj.SegRanges);
        end
    end

    methods
        function set.Baseline(obj, newBaseline)
            obj.BaselineSetMethod = "manual";
            obj.Baseline_I = newBaseline;
        end

        function baseline = get.Baseline(obj)
             baseline = obj.Baseline_I;
        end

        function set.ColorAboveBaseline(obj, color)
            obj.ColorAboveBaselineSetMethod = "manual";
            obj.ColorAboveBaseline_I = validatecolor(color);
        end

        function colorAboveBaseline = get.ColorAboveBaseline(obj)
            colorAboveBaseline = obj.ColorAboveBaseline_I;
        end

        function set.ColorBelowBaseline(obj, color)
            obj.ColorBelowBaselineSetMethod = "manual";
            obj.ColorBelowBaseline_I = validatecolor(color);
        end

        function colorBelowBaseline = get.ColorBelowBaseline(obj)
            colorBelowBaseline = obj.ColorBelowBaseline_I;
        end

    end

end

function [x, y] = validateInputs(x, y, labels)
    x_size = size(x);
    y_size = size(y);

    if x_size(1) == 0 || x_size(2) == 0 || y_size(1) == 0 || y_size(2) == 0
        error("Horizon chart cannot be constructed with empty data");
    end

    if ~isreal(x) || ~isreal(y)
        error("Chart does not work for complex data")
    end 

    if x_size(1) ~= y_size(1)
        error("Number of datapoints for each slice does not match between X and Y");
    end

    if ~validateIsIncreasing(x)
        error("X values should be strictly monotonically increasing")
    end

    if x_size(2) > 1 && x_size(2) < y_size(2)
        error("Number of slices for X-Data can either be 1 or need to match the Y-Data");
    end

    if numel(labels) ~= 0 && numel(labels) ~= y_size(2)
        error("Size of Labels is incorrect. It should either be empty or equal to the number of slices/ 1st dimension of Y data");
    end

end


function isIncreasing = validateIsIncreasing(x)
    for slice = 1:size(x, 2)
        if ~issorted(x(:, slice), 'strictascend')
            isIncreasing = false;
            return;
        end
    end

    isIncreasing = true;
end

function resultAr = binData(data, bins)
    resultAr = histcounts(data, bins);
end

function [x_data, y_data] = transformInputData(x_data, y_data)
    if size(x_data, 1) == 1
        x_data = x_data';
    end

    if size(y_data, 1) == 1
        y_data = y_data';
    end
end 

function [x_vertices, y_vertices] = generatePolygonPoints(dataX, dataY, lower, upper, isSegmentOverBaseline, containsMappingToolbox)

    if isSegmentOverBaseline
        keep = dataY >= lower;
    else
        keep = dataY < upper;
    end
    
    xi = dataX(keep);
    yi = dataY(keep);
    
    yi(yi >= upper) = upper;
    yi(yi < lower) = lower;

    if containsMappingToolbox
        [x_u, y_u] = polyxpoly(dataX, dataY, dataX, upper * ones(1, numel(dataX)));
        [x_l, y_l] = polyxpoly(dataX, dataY, dataX, lower * ones(1, numel(dataX)));
        xi = [xi, x_u', x_l'];
        yi = [yi, y_u', y_l'];
    end


    [xi, idx] = sort(xi);
    yi = yi(idx);

    x_vertices = [xi fliplr(xi)];
    if isSegmentOverBaseline
        y_vertices = [yi ones(1, numel(xi)) * lower];
    else
        y_vertices = [yi ones(1, numel(xi)) * upper];
    end
end

function yi = transformPolygon(yi, lower, upper, baseline)
        
    if lower >= baseline
        yi = yi - (lower - baseline);
    else
        yi = yi + (baseline - upper);
        yi = baseline + abs(baseline - yi);
    end
end

