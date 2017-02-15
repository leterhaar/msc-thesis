function add_14(filename)
output = fopen('.14tmp', 'w');
f = fopen(filename);


while not(feof(f))
    no = str2double(fgetl(f));
    fprintf(output, '%i\n', no+14);
end
fclose(f);
fclose(output);
movefile('.14tmp', filename);
end
