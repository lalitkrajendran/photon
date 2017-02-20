function [num_valid_vectors, num_total_vectors] = count_valid_vectors(data)

    % this is the total number of vectors
    num_total_vectors = size(data.C(:,:,1), 1) * size(data.C(:,:,1), 2);
    
    % calculate the peak ratio
    peak_ratio = data.C(:,:,2) ./ data.C(:,:,3);
    
    % this array contains ones where peak_ratio > 1.2 and zero otherwise
    validity_mask = peak_ratio >= 1.2;
    
    % count the number of non-zero elements
    num_valid_vectors = sum(validity_mask(:));
    

end