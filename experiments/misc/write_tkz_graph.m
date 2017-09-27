function write_tikz_graph(A,Axy,labels,values,directed)

n= size(A);
if isempty(labels), 
    labels= arrayfun(@(x) num2str(x), 1:n, 'UniformOutput', false);
end

if directed
    % bidir
    B= (A&A').*A;
    [Bi Bj Bv] = find(B);
    [Ai Aj Av] = find(A-B);
else
    [Ai Aj Av] = find(triu(A));
end

for i=1:n
    fprintf('\\Vertex[x= %8g, y= %8g]{%s}\n', Axy(i,1), Axy(i,2), labels{i});
end

if directed
    for i=1:length(Ai)
        if values, valstr= sprintf('[label=%g]', Av(i)); else valstr= ''; end
        fprintf('\\Edge%s(%s)(%s)\n', valstr, labels{Ai(i)}, labels{Aj(i)});
    end
    for i=1:length(Bi)
        if values, valstr= sprintf('[label=%g]', Bv(i)); else valstr= ''; end
        fprintf('\\Edge%s(%s)(%s)\n', valstr, labels{Bi(i)}, labels{Bj(i)});
    end
else
    for i=1:length(Ai)
        if values, valstr= sprintf('[label=%g]', Av(i)); else valstr= ''; end
        fprintf('\\Edge%s(%s)(%s)\n', valstr, labels{Ai(i)}, labels{Aj(i)});
    end
end




