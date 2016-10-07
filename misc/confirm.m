function confirm(question)

    key = input([question ' [Y/n]? '], 's');

    if key == 'n'
        error('Not confirmed');

    end 

end