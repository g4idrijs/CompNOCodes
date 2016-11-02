% Save a big list of complementary pairs of different lengths
% (for testing effects of length)

% Code lengths fo generate
lens = 101:300;
for len = lens
    x = genCompPair( ones(1,len-1) );
    fileDir = 'C:\Users\User\Google Drive\Grad_School\HighSpeedCodes\CompNOCodes\Complementary Pairs\PR_CompPairs_VarLen';
    save(strcat(fileDir,'\pairLen',num2str(len),'.mat'),'x');
end
