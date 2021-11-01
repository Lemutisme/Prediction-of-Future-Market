function ezwrite(fname,out_struct,varargin)
% Write a CSV file that contains both text and numeric columns
% from the contents of a structure (of the format produced by "EZread"), 
% without having to manually specify the column type using %f, %s, etc.
%
%   Inputs:
%      1.  Filename
%      2.  structure to write to the file (vector of data per field in
%      structure, as produced by EZread).
%
%   Outputs
%      - none
%
% Paul Taylor
% The MathWorks Australia Pty Ltd
% 20/02/2007

% Set the file delimiter
if nargin == 3
    file_delim = varargin{1};
else
    file_delim = ',';
end

% Open the file
fid = fopen(fname,'w');

if fid == -1
    % could not open the file
    error('Error opening file')
end

headers = fieldnames(out_struct);

% Write header line
for i = 1:length(headers)
    fprintf(fid,headers{i});
    if i < length(headers)
        fprintf(fid,file_delim);
    end
end
fprintf(fid,'\n');

% Determine the data types for each field
format_cell = cell(length(headers),1);
for i = 1:length(headers)
    if isnumeric(out_struct.(headers{i})(1))
        format_cell{i} = 'N';
    else
        format_cell{i} = 'T';
    end
end

% Write the data

for i = 1:length(out_struct.(headers{1}))
    for j = 1:length(headers)
        switch format_cell{j}
            case 'N'
                fprintf(fid,'%g',out_struct.(headers{j})(i));
            case 'T'
                fprintf(fid,'%s',out_struct.(headers{j}){i});
        end

        if j < length(headers)
            fprintf(fid,file_delim);
        end
    end
    fprintf(fid,'\n');
end

% Close the file
fclose(fid);
