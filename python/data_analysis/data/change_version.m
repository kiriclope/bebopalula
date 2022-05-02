function change_version(myDir, name)

  fprintf(myDir)
  fprintf('\n')
  myFiles = dir( fullfile(myDir, '*.mat') );
	
  fprintf(int2str(length(myFiles))) 
  fprintf('\n')

  for k = 1:length(myFiles)
      baseFileName = myFiles(k).name;
      fullFileName = fullfile(myDir, baseFileName);
      fprintf(fullFileName)
      fprintf('\n')
      dum = load(fullFileName);
      newFileName = fullfile(myDir, strcat('new_', baseFileName) );
      % newFileName = fullfile(myDir, strcat(name, '_day_', int2str(k), '.mat') );
      fprintf(newFileName) 
      fprintf('\n')
      save(newFileName, '-struct', 'dum', '-v7');
  end
end
