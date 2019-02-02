function computeripserfiltration(points, x0_index, dx)

    import edu.stanford.math.plex4.*;
    %{
    Here we take the distance matrix and use it to compute the filtration
    values of simplices. In the eccentricity filtration, we only need to 
    compute the filtration values for 0 and 1-simplices as this filtration
    is dependent on only the 1-skeleton. If you change the filtration, you
    may need to explicitly add the filtration values for higher dimensional
    simplices.
    %}
    
    n = size(points);
    n = n(2);
    %the matrix of 1-dimensional filtration values. Note entries (i,i) 
    %represent filtration values of the vertices.
    filtd = zeros(n);
    %find the eccentricity of the seleted basepoint
    ecc = max(dx(:,x0_index));
    %compute filtration values into a matrix
    for i=1:n
        for j=(i):n
            filtd(i,j) = filtration(x0_index, dx, [i j], ecc);
        end
    end
    filtd = filtd + filtd';
    for i=1:n
        filtd(i,i) = filtd(i,i)/2;
    end
    
    stream = api.Plex4.createExplicitSimplexStream(10);
    %add vertices at given time of arrival
    for i=1:n
        stream.addVertex(i, max(min(filtd(i,:)), .5 * (ecc-filtd(i, x0_index))));
    end
    %add 1-simplices at given time of arrival
    for i=2:n
        for j=1:i-1
            stream.addElement([j,i],filtd(j,i));
        end
    end
    %{
    This is an example of what needs to be done if higher dimensional
    simplices time of arrival needs to be added explicitly. Here, this code
    adds the 2 dimensional simplices. Note that having to perform this in 
    Javaplex in dimensions 2 and higher becomes incredibly inneficient, and
    drastically reduces the number of points you can use.
    
    for i=3:n
        for j=2:i-1
            for k=3:j-1
                stream.addElement([k,j,i],filtration(x0_index, dx, [k j i],ecc));
            end
        end
    end
    %}
    stream.finalizeStream();
    
    %Set the maximum dimension to compute barcodes using Javaplex. 
    %Setting max_dimension=n will cause computation of barcodes in up
    %to dimension n-1.
    max_dimension = 1;
    
    %sets the maximum filtration value for Javaplex to use to compute.
    %value should be changed based on size/diameter of particular dataset
    %to have optimal visualization of barcodes
    max_filtration_value = 70;


    % get persistence algorithm over Z/2Z, can substitute 2 for any prime p
    persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);

    % compute the intervals
    intervals = persistence.computeIntervals(stream);

    % create the barcode plots
    options.filename = 'Local Filtration';
    options.max_filtration_value = max_filtration_value;

    plot_barcodes(intervals, options);
    
    %Check to determine if code is running/stalling:
    fprintf('we have reached ripser computation\n');
    
    %Use ripser to compute homology in higher dimensions
    Is2geo = ripserDM(filtd, 2, 2);
    disp(Is2geo{2});
    
    %plot results
    figure

    plotDGM(Is2geo{2}, 'r', 20, 0), axis square;
    hold on;
    plot(0:.01:50, 0:0.01:50, 'k');
    %axis limits should be changed to give best viewing window for 
    %a particular dataset, varies greatly depending on diameter of 
    %dataset in R^3.
    xlim([0,20])
    ylim([0,20])
    
    

end


%{
This is the function which computes filtration values. It is currently
coded to compute this for the eccentricity filtration, but can be adapted
to compute the filtration values for any local basepoint filtration.
%}

function filtration_value = filtration(x0_index, dW, indices, ecc)
    pairdist = [];
    for j=1:length(indices)
        for k=(j+1):length(indices)
           pairdist = [pairdist; dW(indices(k),indices(j))];
        end
    end
    
    basedist = [];
    for j=1:length(indices)
        basedist = [basedist; dW(x0_index, indices(j))];
    end
    
    filtration_value = max(max(pairdist), 0.5*(ecc - min(basedist)));
    
end
