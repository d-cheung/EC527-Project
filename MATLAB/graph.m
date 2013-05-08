function graph( input_array, label, color )
% Graphs the input array with the label label and the color color
%   The input array consists of the following per image:
%       Image Width, Image Height, Detected Features, Cumulative Mean
%       Image Width, Image Height, Detected Features, IntegralImage Mean
%       Image Width, Image Height, Detected Features, FastHessian Mean
%       Image Width, Image Height, Detected Features, SurfDescriptor Mean
%
%   Figure 1 contains the Total Time vs. Image Size
%   Figure 2 contains the IntegralImage Time vs. Image Size
%   Figure 3 contains the FastHessian Time vs. Image Size
%   Figure 4 contains the Surf Descriptor Time vs. Number of Interest Points

    for i=0:((size(input_array,1) / 4) - 1)
             totalSize(i+1) = input_array(4*i+1,3);
         integralImage(i+1) = input_array(4*i+2,3);
           fastHessian(i+1) = input_array(4*i+3,3);
        surfDescriptor(i+1) = input_array(4*i+4,3);
              elements(i+1) = input_array(4*i+1,1);
        interestPoints(i+1) = input_array(4*i+1,2);
    end

    
    figure(1);
    hold on;
    [Y, I] = sort(elements);
    scatter(elements, totalSize, 100, color, 'DisplayName', label);
    plot(elements(I), totalSize(I), '-', 'Color', color, 'DisplayName', '');
    xlabel('Image Size (Pixels)', 'FontSize', 20);
    ylabel('Time (ms)', 'FontSize', 20);
    title('Total Time vs. Image Size', 'FontSize', 20);
    
    figure(2);
    hold on;
    scatter(elements, integralImage, 100, color, 'DisplayName', label);
    plot(elements(I), integralImage(I), '-', 'Color', color, 'DisplayName', '');
    xlabel('Image Size (Pixels)', 'FontSize', 20);
    ylabel('Time (ms)', 'FontSize', 20);
    title('IntegralImage Time vs. Image Size', 'FontSize', 20);
    
    figure(3);
    hold on;
    scatter(elements, fastHessian, 100, color, 'DisplayName', label);
    plot(elements(I), fastHessian(I), '-', 'Color', color, 'DisplayName', '');
    xlabel('Image Size (Pixels)', 'FontSize', 20);
    ylabel('Time (ms)', 'FontSize', 20);
    title('FastHessian Time vs. Image Size', 'FontSize', 20);

    figure(4);
    hold on;
    [Y, I] = sort(interestPoints);
    scatter(interestPoints, surfDescriptor, 100, color, 'DisplayName', label);
    plot(interestPoints(I), surfDescriptor(I), '-', 'Color', color, 'DisplayName', '');
    xlabel('Number of Interest Points', 'FontSize', 20);
    ylabel('Time (ms)', 'FontSize', 20);
    title('Surf Descriptor Time vs. Number of Interest Points', 'FontSize', 20);

end

