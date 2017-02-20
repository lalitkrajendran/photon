%% =========================================================
% Get the meta data as numbers from nrrd
%  It is useful to use with:
%  http://de.mathworks.com/matlabcentral/fileexchange/34653-nrrd-format-file-reader
%  http://de.mathworks.com/matlabcentral/fileexchange/48621-nrrdwriter-filename--matrix--pixelspacing--origin--encoding-
%
% =========================================================
 function y=nrrdmeta2num(imgh)
% problem is all the fields are strings so we need to extract numbers
% there should be simpler way using regular expression ;)
% Example of nrrd header
% imgh =         type: 'short'
%           dimension: '3'
%               space: 'left-posterior-superior'
%               sizes: '185 109 53'
%     spacedirections: '(0.125,0,0) (0,0.125,0) (0,0,0.375)'
%               kinds: 'domain domain domain'
%              endian: 'little'
%            encoding: 'gzip'
%         spaceorigin: '(-4.1875,-9.6875,-10.125)'
% 
% to use nrrdwriter we need:  pixelspacing, origin, encoding
%          pixelspacing: extracted from spacedirections   
%                origin: spaceorigin
%              encoding: encoding
%
%   nrrdwriter(filename, 3d image, pixelspacing, origin, encoding)
%   example: 
%   nrrdWriter('mynrrd.nrrd',x,[0.125,0.125,0.375],[-4.1875,-9.6875,-10.125],'gzip');
% example for use: 
%  [im imgh]=nrrdread('My_3D_Image.nrrd'); 
%   y=nrrdmeta2num(imgh);
%   nrrdwriter('My_3D_Image_Copy.nrrd',im,y{1},y{2},y{3});
spcs=imgh.spacedirections;
orgs=imgh.spaceorigin;
encs=imgh.encoding;
  
p1 = findstr('(', spcs); p2 = findstr(')', spcs); Nms= findstr(',', spcs);
s1=str2num(spcs(p1(1)+1:Nms(1)-1)) ; 
s2=str2num(spcs(Nms(3)+1:Nms(4)-1));
s3=str2num(spcs(Nms(6)+1:p2(3)-1)) ;
spc=[s1,s2,s3];

p1 = findstr('(', orgs); p2 = findstr(')', orgs); Nms= findstr(',', orgs);
s1=str2num(orgs(p1(1)+1:Nms(1)-1)) ;
s2=str2num(orgs(Nms(1)+1:Nms(2)-1));
s3=str2num(orgs(Nms(2)+1:p2(1)-1)) ;
org=[s1,s2,s3];


y= {spc org encs};
