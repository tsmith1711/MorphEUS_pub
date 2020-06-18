function [outputArg1] = printstructfun(x,transverse_file)
%function previously contained in MicrobeJ_segmentation with the addition
%of a timing function that gives the time taken to process each drug at
%what is the longest section of the code 

 outputArg1 = structfun(@(y) ... %for each photo in bacteria
                        cellfun(@(z) ... %for each NAME_id in bacteria
                            transverse_file(strcmp(transverse_file.NAME_id, char(z)), :), ... %find rows of matching id in transverse
                        y.NAME_id, 'UniformOutput', false).', ...
                 x, 'UniformOutput', false);

fprintf('Processed a drug. Total time elapsed: %.2fs \n',toc)
end

