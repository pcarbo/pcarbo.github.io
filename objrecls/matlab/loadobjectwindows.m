% LOADOBJECTWINDOWS loads the object windows from a file.
%   LOADOBJECTWINDOWS(F) returns an array of object windows by reading
%   the data from a file specified by pathname F.

function objs = loadobjectwindows (f)

  % Open the file.
  infile = fopen(f, 'r');
  if infile == -1,
    error(sprintf('Unable to open file %s for reading', f));
  end;

  objs = [];  % The object structure.
  s    = fgetl(infile);
  while ischar(s),
    [t s] = strtok(s);
    [t s] = strtok(s);
    obj.x = max(0,str2num(t));
    [t s] = strtok(s);
    obj.y = max(0,str2num(t));
    [t s] = strtok(s);
    obj.w = str2num(t);
    obj.h = str2num(s);
    s     = fgetl(infile);
    
    objs = [objs obj];
  end;

  % Close the file.
  fclose(infile);    
