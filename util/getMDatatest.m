function data = getMDatatest ( n )
ImageSize = 256; 
TesDataNumber = 50 ;

data.train = single(zeros(ImageSize, ImageSize));
data.label = single (zeros(ImageSize, ImageSize));
%%sampling rate
dir = '.\data_20\test\';
tep = load (strcat(dir , saveName(n,floor(log10(TesDataNumber))+1)));
data.label = tep.data.label;
data.train = tep.data.train;