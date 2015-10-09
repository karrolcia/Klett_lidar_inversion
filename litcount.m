% The function litcount determines whether an input literal string (literal) 
% appears in each line. If it does, the function prints the entire line preceded by 
% the number of times the literal string appears on the line. If it does, the function 
% prints the entire line preceded by the number of times the literal string appears on the line.


function y = litcount(filename, literal)

fid = fopen(filename);
y = 0;
tline = fgetl(fid);
while ischar(tline)
   matches = strfind(tline, literal);
   num = length(matches);
   if num > 0
      y = y + num;
      fprintf(1,'%d:%s\n',num,tline);
   end
   tline = fgetl(fid);
end
fclose(fid);