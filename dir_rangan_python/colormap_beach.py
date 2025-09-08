import numpy as np

'''
function output = colormap_beach(n_cra);
if (nargin<1); n_cra = 64; end;
if n_cra<0; imagesc(1:64); colormap(colormap_beach()); return; end;
gamma1 = 0.5;
gamma2 = 0.25;
cra = zeros(n_cra,3);
m1 = floor(n_cra/2);
cra(m1,:) = [1,1,1];
tmp_ = linspace(0,1,m1-1).^gamma1;
cra(1:m1-1,1) = transpose(tmp_).^6;
cra(1:m1-1,2) = 0.5 + 0.5*transpose(tmp_);
cra(1:m1-1,3) = ones(m1-1,1);
m2 = ceil(n_cra/2);
cra(m2,:) = [1,1,1];
tmp_ = linspace(1,0,n_cra-m2).^gamma2;
cra(m2+1:end,1) = ones(n_cra-m2,1);
cra(m2+1:end,2) = transpose(tmp_).^6;
cra(m2+1:end,3) = 0.5 + 0.5*transpose(tmp_);
cra(:,3) = cra(:,3).*transpose(sqrt(sqrt(linspace(-1,+1,n_cra).^2)));
output = cra;
'''
def colormap_beach(n_cra=64):
    if n_cra < 0:
        import matplotlib.pyplot as plt
        plt.imshow(np.arange(64).reshape(1, -1), aspect='auto', cmap=plt.cm.get_cmap('viridis'))
        plt.colorbar()
        plt.show()
        return

    gamma1 = 0.5
    gamma2 = 0.25
    cra = np.zeros((n_cra, 3))
    m1 = n_cra // 2
    cra[m1, :] = [1, 1, 1]
    tmp_ = np.linspace(0, 1, m1 - 1) ** gamma1
    cra[:m1 - 1, 0] = tmp_ ** 6
    cra[:m1 - 1, 1] = 0.5 + 0.5 * tmp_
    cra[:m1 - 1, 2] = 1
    m2 = (n_cra + 1) // 2
    cra[m2 - 1, :] = [1, 1, 1]
    tmp_ = np.linspace(1, 0, n_cra - m2) ** gamma2
    cra[m2:, 0] = 1
    cra[m2:, 1] = tmp_ ** 6
    cra[m2:, 2] = 0.5 + 0.5 * tmp_
    cra[:, 2] *= np.sqrt(np.sqrt(np.linspace(-1, 1, n_cra) ** 2))
    return cra
