function output_data = reorder_threads( input_data )

    t = 1;
    for i=1:5:3*5
        output_data(t+0,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+1,:))];
        output_data(t+1,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+2,:))];
        output_data(t+2,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+3,:))];
        output_data(t+3,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+4,:))];
        t = t + 4;
    end
    for i=16:5:size(input_data,1)
        output_data(t+0,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+1,1:10))];
        output_data(t+1,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+2,1:10))];
        output_data(t+2,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+3,1:10))];
        output_data(t+3,:) = [input_data(i,1) * input_data(i,2), input_data(i,3), input_data(i, 4), mean(input_data(i+4,1:10))];
        t = t + 4;
    end
        
        
end

