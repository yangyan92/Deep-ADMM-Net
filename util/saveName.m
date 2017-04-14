function name = saveName(n, l)
name = num2str(n);
for i=1:l-length(name)
    name = ['0' name];
end;

end

