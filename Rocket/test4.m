D = parallel.pool.DataQueue;
h = waitbar(0, 'Please wait ...');
afterEach(D, @nUpdateWaitbar);

N = 200;
p = 1;

parfor (i = 1:N, 0)
    a(i) = max(abs(eig(rand(400))));
    send(D, i);
end

function nUpdateWaitbar(~)
    waitbar(p/N, h);
    p = p + 1;
end