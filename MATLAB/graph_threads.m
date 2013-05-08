function graph_threads( input_array )



    for j=0:5
        for i=0:7
                 totalSize(i+1) = input_array(j*32+4*i+1,4);
             integralImage(i+1) = input_array(j*32+4*i+2,4);
               fastHessian(i+1) = input_array(j*32+4*i+3,4);
            surfDescriptor(i+1) = input_array(j*32+4*i+4,4);
                  elements(i+1) = input_array(j*32+4*i+1,1);
            interestPoints(i+1) = input_array(j*32+4*i+1,3);
        end
        
        threads = input_array(j*32+1, 2);
        
        figure(5);
        hold on;
        scatter(elements, totalSize, 'DisplayName', num2str(threads));

        figure(6);
        hold on;
        scatter(elements, integralImage, 'DisplayName', num2str(threads));

        figure(7);
        hold on;
        scatter(elements, fastHessian, 'DisplayName', num2str(threads));

        figure(8);
        hold on;
        scatter(interestPoints, surfDescriptor, 'DisplayName', num2str(threads));
    end

end

