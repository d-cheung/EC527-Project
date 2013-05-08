function output_data = reorder( input_data )
% Reorders array for graphability
%   One image in the the input data is in the following format:
%       Image Width, Image Height, Detected Features
%       Cumulative times: t1, t2 ... tx
%       Integral Image times: t1, t2 ... tx
%       FastHessian times: t1, t2 ... tx
%       SurfDescriptor times: t1, t2 ... tx
%   This is repeated for every image thereafter.
%
%   The output array consists of the following per image:
%       Image Width, Image Height, Detected Features, Cumulative Mean
%       Image Width, Image Height, Detected Features, IntegralImage Mean
%       Image Width, Image Height, Detected Features, FastHessian Mean
%       Image Width, Image Height, Detected Features, SurfDescriptor Mean

    t = 1;
    for i=1:5:size(input_data,1)
        output_data(t+0,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), mean(input_data(i+1,:))];
        output_data(t+1,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), mean(input_data(i+2,:))];
        output_data(t+2,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), mean(input_data(i+3,:))];
        output_data(t+3,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), mean(input_data(i+4,:))];
        t = t + 4;
    end
        
        
end

