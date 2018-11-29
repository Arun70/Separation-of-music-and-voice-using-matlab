[x, fs] = wavread('105 - The Battle Over the Barrier_001 (online-audio-converter.com).wav');
len = 0.04;                            
N = 2 .^ nextpow2(len * fs);			
window = hamming(N, 'periodic');		
step = N / 2;                          
[t, k] = size(x);                      
if (k == 2)
    x = x(:, 1) + x(:, 2);             
end
k = 1;
X = stft(x(:, 1), window, step); 
V = abs(X(1 : N / 2 + 1, 1, 1));        
VSquare = mean(V .^ 2, 3);
b = beatSpectrum(V);
b(1) = [];                              
up = min(8, (length(x) / fs) / 3);
range = [0.8, up];
range = ceil(range * fs / step);
b = b(range(1) : range(2));
[~, p] = max(b);                        
p = p + range(1) - 1;
y = zeros(t, k);                        
vocal = zeros(t, k); 
Mask = mask(V(:, :, 1), p);
Mask = cat(1, Mask, flipud(Mask(1 : end - 2, :)));
yPart = istft(Mask .* X(:, :, 1), window, step);
y(:, 1) = yPart(1 : t);
vocal = x - y; 