% Load the data from the text file
data = readmatrix('k-mer_sequence.txt');

% Initialize an empty matrix to store the results
results = [];

% Use a for loop to process the data
for i = 1:length(data)

    [RNAbracket,Energy] = rnafold(i);

    % Append the result to the results matrix
    results = [results; Energy];
end

% Write the results to a new text file
writematrix(results, 'k-mer_MFE-output.txt');