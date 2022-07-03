function [stop, options, optchanged] =myoutSA(options, optimvalues, flag)
    global fitall_sa
    optchanged=false;
    stop=false;
    switch flag
        case 'init'
            disp('starting Simulated annealing algorithm...');
        case 'iter'
            fitall_sa(optimvalues.iteration)=optimvalues.bestfval;
            fprintf('\nIteration #%d \t error = ', optimvalues.iteration);
            disp(optimvalues.bestfval);
            fprintf('\t\nx = ');
            disp(optimvalues.x);
            fprintf('\t\nbestx = ');
            disp(optimvalues.bestx);
        case 'done'
            disp('ending SA...');
    end
end