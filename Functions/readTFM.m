function mat = readTFM(tfm)
% function to read tfm

fid = fopen(tfm);



for i = 1:4
    %
    newdat = [];
    
    %
    ln = fgets(fid);
    if contains(ln, sprintf('\t'))
        data = strip(split(ln,sprintf('\t')));
    elseif contains(ln, sprintf('  '))
        data = strip(split(ln,'  '));
    else
        data = strip(split(ln,' '));
    end
    for j = 1:length(data)
        if ~isempty(data{j})
            newdat(1,end+1) = str2double(data{j});
        end
    end
    mat(i,:) = newdat';
end
fclose('all');