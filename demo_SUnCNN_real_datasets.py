# -*- coding: utf-8 -*-
"""
用于测试ADSU算法  Cuprite
"""
from __future__ import print_function
import matplotlib.pyplot as plt
#%matplotlib inline

import os
#os.environ['CUDA_VISIBLE_DEVICES'] = '3'
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import numpy as np
from models import *
import sys
import torch
import torch.optim

from skimage.metrics import peak_signal_noise_ratio as compare_psnr
from skimage.metrics import mean_squared_error as compare_mse
from utils.denoising_utils import *

from skimage._shared import *
from skimage.util import *
from skimage.metrics.simple_metrics import _as_floats
from skimage.metrics.simple_metrics import mean_squared_error


from UtilityMine import *
from sr_utils import tv_loss
from numpy import linalg as LA

torch.backends.cudnn.enabled = True
torch.backends.cudnn.benchmark =True
dtype = torch.cuda.FloatTensor

PLOT = False
import scipy.io

#%% 替换数据即可
fname2  = "ADSU_Data/Cuprite/image.mat"
mat2 = scipy.io.loadmat(fname2)
img_np_gt = mat2["image"]
img_np_gt = img_np_gt.transpose(2,0,1)
[p1, nr1, nc1] = img_np_gt.shape

# # 加载参考丰度
# fname3  = "ADSU_Data/DC1/XT.mat"
# mat3 = scipy.io.loadmat(fname3)
# A_true_np = mat3["XT"]
# A_true_np = A_true_np.transpose(2,0,1)

# 加载光谱库
fname4  = "ADSU_Data/Cuprite/M.mat"
mat4 = scipy.io.loadmat(fname4)
EE = mat4["M"]

#%% 光谱库中端元数
LibS=EE.shape[1]
save_result = True

import time
from tqdm import tqdm
# rmax = 5;
tol1 = 1

for fi in tqdm(range(tol1)):
    img_noisy_np = img_np_gt
    INPUT = 'noise'  # 'meshgrid'
    pad = 'reflection'
    need_bias = True
    OPT_OVER = 'net'  # 'net,input'

    reg_noise_std = 0.0
    LR1 = 0.001

    OPTIMIZER1 = 'adam'  # 'RMSprop'#'adam' # 'LBFGS'
    show_every = 500
    exp_weight = 0.99

    if fi == 0:
        num_iter1 = 20000
    input_depth =img_noisy_np.shape[0]
    class CAE_AbEst(nn.Module):
        def __init__(self):
            super(CAE_AbEst, self).__init__()
            self.conv1 = nn.Sequential(
                UnmixArch(
                    input_depth, EE.shape[1],
                    num_channels_down=[256],
                    num_channels_up=[256],
                    num_channels_skip=[4],
                    filter_size_up=3, filter_size_down=3, filter_skip_size=1,
                    upsample_mode='bilinear',  # downsample_mode='avg',
                    need1x1_up=True,
                    need_sigmoid=True, need_bias=True, pad=pad, act_fun='LeakyReLU').type(dtype)
            )

        def forward(self, x):
            x = self.conv1(x)
            return x


    net1 = CAE_AbEst()
    net1.cuda()
    print(net1)

    # Compute number of parameters
    s = sum([np.prod(list(p11.size())) for p11 in net1.parameters()]);
    print('Number of params: %d' % s)

    # Loss
    mse = torch.nn.MSELoss().type(dtype)
    img_noisy_torch = np_to_torch(img_noisy_np).type(dtype)
    # if fk==0:
    net_input1 = get_noise(input_depth, INPUT,
                           (img_noisy_np.shape[1], img_noisy_np.shape[2])).type(dtype).detach()
    net_input1 = img_noisy_torch
    E_torch = np_to_torch(EE).type(dtype)
    # %%
    net_input_saved = net_input1.detach().clone()
    noise = net_input1.detach().clone()
    out_avg = None
    out_HR_avg = None
    last_net = None
    RMSE_LR_last = 0
    loss = np.zeros((num_iter1, 1))
    AE = np.zeros((num_iter1, 1))
    i = 0


    def closure1():

        global i, RMSE_LR, RMSE_LR_ave, RMSE_HR, out_LR_np, out_avg_np, out_LR \
            , out_avg, out_HR_np, out_HR_avg, out_HR_avg_np, RMSE_LR_last, last_net \
            , net_input, RMSE_LR_avg, RMSE_HR_avg, RE_HR_avg, RE_HR, Eest, loss, AE \
            , MAE_LR, MAE_LR_avg, MAE_HR, MAE_HR_avg

        if reg_noise_std > 0:
            net_input = net_input_saved + (noise.normal_() * reg_noise_std)

        out_LR = net1(net_input1)
        out_HR = torch.mm(E_torch.view(p1, LibS), out_LR.view(LibS, nr1 * nc1))
        # Smoothing
        if out_avg is None:
            out_avg = out_LR.detach()
            out_HR_avg = out_HR.detach()
        else:
            out_avg = out_avg * exp_weight + out_LR.detach() * (1 - exp_weight)
            out_HR_avg = out_HR_avg * exp_weight + out_HR.detach() * (1 - exp_weight)

        # %%
        out_HR = out_HR.view((1, p1, nr1, nc1))
        total_loss = mse(img_noisy_torch, out_HR)
        total_loss.backward()
        if 1:
            out_LR_np = out_LR.detach().cpu().squeeze().numpy()
            out_avg_np = out_avg.detach().cpu().squeeze().numpy()

            # SRE = 10 * np.log10(
            #     LA.norm(A_true_np.astype(np.float32).reshape((EE.shape[1], nr1 * nc1)), 'fro') / LA.norm(
            #         (A_true_np.astype(np.float32) - np.clip(out_LR_np, 0, 1)).reshape((EE.shape[1], nr1 * nc1)), 'fro'))
            # SRE_avg = 10 * np.log10(
            #     LA.norm(A_true_np.astype(np.float32).reshape((EE.shape[1], nr1 * nc1)), 'fro') / LA.norm(
            #         (A_true_np.astype(np.float32) - np.clip(out_avg_np, 0, 1)).reshape((EE.shape[1], nr1 * nc1)),
            #         'fro'))
            # MAE_LR = 100 * np.mean(abs(A_true_np.astype(np.float32) - np.clip(out_LR_np, 0, 1)))
            # MAE_LR_avg = 100 * np.mean(abs(A_true_np.astype(np.float32) - np.clip(out_avg_np, 0, 1)))
            # print('Iteration %05d    Loss %f   MAE_LR: %f MAE_LR_avg: %f  SRE: %f SRE_avg: %f' % (
            # i, total_loss.item(), MAE_LR, MAE_LR_avg, SRE, SRE_avg), '\r', end='')

            print('Iteration %05d    Loss %f' % (i, total_loss.item()), '\r', end='')
        if PLOT and i % show_every == 0:
            out_LR_np = torch_to_np(out_LR)
            out_avg_np = torch_to_np(out_avg)
            plt.plot(out_LR_np.reshape(LibS, nr1 * nc1))
        loss[i] = total_loss.item()
        i += 1

        return total_loss


    p11 = get_params(OPT_OVER, net1, net_input1)
    optimize(OPTIMIZER1, p11, closure1, LR1, num_iter1)
    if 1:
        out_LR_np = out_LR.detach().cpu().squeeze().numpy()
        out_avg_np = out_avg.detach().cpu().squeeze().numpy()

        # MAE_LR_avg = 100 * np.mean(abs(A_true_np.astype(np.float32) - np.clip(out_avg_np, 0, 1)))
        # MAE_LR = 100 * np.mean(abs(A_true_np.astype(np.float32) - np.clip(out_LR_np, 0, 1)))
        # SRE = 10 * np.log10(LA.norm(A_true_np.astype(np.float32).reshape((EE.shape[1], nr1 * nc1)), 'fro') / LA.norm(
        #     (A_true_np.astype(np.float32) - np.clip(out_LR_np, 0, 1)).reshape((EE.shape[1], nr1 * nc1)), 'fro'))
        # SRE_avg = 10 * np.log10(
        #     LA.norm(A_true_np.astype(np.float32).reshape((EE.shape[1], nr1 * nc1)), 'fro') / LA.norm(
        #         (A_true_np.astype(np.float32) - np.clip(out_avg_np, 0, 1)).reshape((EE.shape[1], nr1 * nc1)), 'fro'))
        # print('Iteration %05d  MAE_LR: %f MAE_LR_avg: %f  SRE: %f SRE_avg: %f ' % (i, MAE_LR, MAE_LR_avg, SRE, SRE_avg),
        #       '\r', end='')

        print('Iteration %05d' % (i),
              '\r', end='')
    if  save_result is True:
        scipy.io.savemat(r"D:\work\TGRS\data\Cuprite\out_avg_np%01d%01d.mat" % (fi+2,1),
                         {'out_avg_np%01d%01d' % (fi+2,1):out_avg_np.transpose(1,2,0)})
        scipy.io.savemat(r"D:\work\TGRS\data\Cuprite\out_LR_np%01d%01d.mat" % (fi+2,1),
                         {'out_LR_np%01d%01d' % (fi+2, 1):out_LR_np.transpose(1,2,0)})


# -*- coding: utf-8 -*-
"""
用于测试ADSU算法 Urban_R162数据

"""

from __future__ import print_function
import matplotlib.pyplot as plt
#%matplotlib inline

import os
#os.environ['CUDA_VISIBLE_DEVICES'] = '3'
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import numpy as np
from models import *
import sys
import torch
import torch.optim

from skimage.metrics import peak_signal_noise_ratio as compare_psnr
from skimage.metrics import mean_squared_error as compare_mse
from utils.denoising_utils import *

from skimage._shared import *
from skimage.util import *
from skimage.metrics.simple_metrics import _as_floats
from skimage.metrics.simple_metrics import mean_squared_error


from UtilityMine import *
from sr_utils import tv_loss
from numpy import linalg as LA

torch.backends.cudnn.enabled = True
torch.backends.cudnn.benchmark =True
dtype = torch.cuda.FloatTensor

PLOT = False
import scipy.io
#%% 加载数据 方块75*75数据(图像数据)
fname2  = "ADSU_Data/real_datasets/Urban_R162.mat"
mat2 = scipy.io.loadmat(fname2)
img_np_gt = mat2["Y"]
img_np_gt = img_np_gt.transpose(2,0,1)
[p1, nr1, nc1] = img_np_gt.shape

# 加载光谱库
fname4  = "ADSU_Data/real_datasets/EE_Urban.mat"
mat4 = scipy.io.loadmat(fname4)
EE = mat4["EE"]

#%% 光谱库中端元数
LibS=EE.shape[1]
save_result = True

import time
from tqdm import tqdm
# rmax = 5;
tol1 = 1

for fi in tqdm(range(tol1)):
    img_noisy_np = img_np_gt
    INPUT = 'noise'  # 'meshgrid'
    pad = 'reflection'
    need_bias = True
    OPT_OVER = 'net'  # 'net,input'

    #
    reg_noise_std = 0.0
    LR1 = 0.001

    OPTIMIZER1 = 'adam'  # 'RMSprop'#'adam' # 'LBFGS'
    show_every = 500
    exp_weight = 0.99

    if fi == 0:
        num_iter1 = 20000
    input_depth =img_noisy_np.shape[0]
    class CAE_AbEst(nn.Module):
        def __init__(self):
            super(CAE_AbEst, self).__init__()
            self.conv1 = nn.Sequential(
                UnmixArch(
                    input_depth, EE.shape[1],
                    num_channels_down=[256],
                    num_channels_up=[256],
                    num_channels_skip=[4],
                    filter_size_up=3, filter_size_down=3, filter_skip_size=1,
                    upsample_mode='bilinear',  # downsample_mode='avg',
                    need1x1_up=True,
                    need_sigmoid=True, need_bias=True, pad=pad, act_fun='LeakyReLU').type(dtype)
            )

        def forward(self, x):
            x = self.conv1(x)
            return x


    net1 = CAE_AbEst()
    net1.cuda()
    print(net1)

    # Compute number of parameters
    s = sum([np.prod(list(p11.size())) for p11 in net1.parameters()]);
    print('Number of params: %d' % s)

    # Loss
    mse = torch.nn.MSELoss().type(dtype)
    img_noisy_torch = np_to_torch(img_noisy_np).type(dtype)
    # if fk==0:
    net_input1 = get_noise(input_depth, INPUT,
                           (img_noisy_np.shape[1], img_noisy_np.shape[2])).type(dtype).detach()
    net_input1 = img_noisy_torch
    E_torch = np_to_torch(EE).type(dtype)
    # %%
    net_input_saved = net_input1.detach().clone()
    noise = net_input1.detach().clone()
    out_avg = None
    out_HR_avg = None
    last_net = None
    RMSE_LR_last = 0
    loss = np.zeros((num_iter1, 1))
    AE = np.zeros((num_iter1, 1))
    i = 0


    def closure1():

        global i, RMSE_LR, RMSE_LR_ave, RMSE_HR, out_LR_np, out_avg_np, out_LR \
            , out_avg, out_HR_np, out_HR_avg, out_HR_avg_np, RMSE_LR_last, last_net \
            , net_input, RMSE_LR_avg, RMSE_HR_avg, RE_HR_avg, RE_HR, Eest, loss, AE \
            , MAE_LR, MAE_LR_avg, MAE_HR, MAE_HR_avg

        if reg_noise_std > 0:
            net_input = net_input_saved + (noise.normal_() * reg_noise_std)

        out_LR = net1(net_input1)
        out_HR = torch.mm(E_torch.view(p1, LibS), out_LR.view(LibS, nr1 * nc1))
        # Smoothing
        if out_avg is None:
            out_avg = out_LR.detach()
            out_HR_avg = out_HR.detach()
        else:
            out_avg = out_avg * exp_weight + out_LR.detach() * (1 - exp_weight)
            out_HR_avg = out_HR_avg * exp_weight + out_HR.detach() * (1 - exp_weight)

        # %%
        out_HR = out_HR.view((1, p1, nr1, nc1))
        total_loss = mse(img_noisy_torch, out_HR)
        total_loss.backward()
        if 1:
            out_LR_np = out_LR.detach().cpu().squeeze().numpy()
            out_avg_np = out_avg.detach().cpu().squeeze().numpy()

            print('Iteration %05d    Loss %f' % (i, total_loss.item()), '\r', end='')
        if PLOT and i % show_every == 0:
            out_LR_np = torch_to_np(out_LR)
            out_avg_np = torch_to_np(out_avg)
            plt.plot(out_LR_np.reshape(LibS, nr1 * nc1))
        loss[i] = total_loss.item()
        i += 1

        return total_loss


    p11 = get_params(OPT_OVER, net1, net_input1)
    optimize(OPTIMIZER1, p11, closure1, LR1, num_iter1)
    if 1:
        out_LR_np = out_LR.detach().cpu().squeeze().numpy()
        out_avg_np = out_avg.detach().cpu().squeeze().numpy()

        print('Iteration %05d' % (i),
              '\r', end='')
    if  save_result is True:
         scipy.io.savemat(r"D:\work\TGRS\data\Urban_R162\Urban\out_avg_np_Urban_%01d%01d.mat" % (fi+2,1),
                        {'out_avg_np_Urban_%01d%01d' % (fi+2,1):out_avg_np.transpose(1,2,0)})
         scipy.io.savemat(r"D:\work\TGRS\data\Urban_R162\Urban\out_LR_np_Urban_%01d%01d.mat" % (fi+2,1),
                        {'out_LR_np_Urban_%01d%01d' % (fi+2, 1):out_LR_np.transpose(1,2,0)})
# %%




# -*- coding: utf-8 -*-
"""
用于测试ADSU算法 jasperRidge2数据

"""

torch.backends.cudnn.enabled = True
torch.backends.cudnn.benchmark =True
dtype = torch.cuda.FloatTensor

PLOT = False
import scipy.io
#%% 加载数据 方块75*75数据(图像数据)
fname2  = "ADSU_Data/real_datasets/jasperRidge2_R198.mat"
mat2 = scipy.io.loadmat(fname2)
img_np_gt = mat2["Y"]
img_np_gt = img_np_gt.transpose(2,0,1)
[p1, nr1, nc1] = img_np_gt.shape

# 加载光谱库
fname4  = "ADSU_Data/real_datasets/EE_jasperRidge2.mat"
mat4 = scipy.io.loadmat(fname4)
EE = mat4["EE"]

#%% 光谱库中端元数
LibS=EE.shape[1]
save_result = True

import time
from tqdm import tqdm
# rmax = 5;
tol1 = 1

for fi in tqdm(range(tol1)):
    img_noisy_np = img_np_gt
    INPUT = 'noise'  # 'meshgrid'
    pad = 'reflection'
    need_bias = True
    OPT_OVER = 'net'  # 'net,input'

    #
    reg_noise_std = 0.0
    LR1 = 0.001

    OPTIMIZER1 = 'adam'  # 'RMSprop'#'adam' # 'LBFGS'
    show_every = 500
    exp_weight = 0.99

    if fi == 0:
        num_iter1 = 20000
    input_depth =img_noisy_np.shape[0]
    class CAE_AbEst(nn.Module):
        def __init__(self):
            super(CAE_AbEst, self).__init__()
            self.conv1 = nn.Sequential(
                UnmixArch(
                    input_depth, EE.shape[1],
                    num_channels_down=[256],
                    num_channels_up=[256],
                    num_channels_skip=[4],
                    filter_size_up=3, filter_size_down=3, filter_skip_size=1,
                    upsample_mode='bilinear',  # downsample_mode='avg',
                    need1x1_up=True,
                    need_sigmoid=True, need_bias=True, pad=pad, act_fun='LeakyReLU').type(dtype)
            )

        def forward(self, x):
            x = self.conv1(x)
            return x


    net1 = CAE_AbEst()
    net1.cuda()
    print(net1)

    # Compute number of parameters
    s = sum([np.prod(list(p11.size())) for p11 in net1.parameters()]);
    print('Number of params: %d' % s)

    # Loss
    mse = torch.nn.MSELoss().type(dtype)
    img_noisy_torch = np_to_torch(img_noisy_np).type(dtype)
    # if fk==0:
    net_input1 = get_noise(input_depth, INPUT,
                           (img_noisy_np.shape[1], img_noisy_np.shape[2])).type(dtype).detach()
    net_input1 = img_noisy_torch
    E_torch = np_to_torch(EE).type(dtype)
    # %%
    net_input_saved = net_input1.detach().clone()
    noise = net_input1.detach().clone()
    out_avg = None
    out_HR_avg = None
    last_net = None
    RMSE_LR_last = 0
    loss = np.zeros((num_iter1, 1))
    AE = np.zeros((num_iter1, 1))
    i = 0


    def closure1():

        global i, RMSE_LR, RMSE_LR_ave, RMSE_HR, out_LR_np, out_avg_np, out_LR \
            , out_avg, out_HR_np, out_HR_avg, out_HR_avg_np, RMSE_LR_last, last_net \
            , net_input, RMSE_LR_avg, RMSE_HR_avg, RE_HR_avg, RE_HR, Eest, loss, AE \
            , MAE_LR, MAE_LR_avg, MAE_HR, MAE_HR_avg

        if reg_noise_std > 0:
            net_input = net_input_saved + (noise.normal_() * reg_noise_std)

        out_LR = net1(net_input1)
        out_HR = torch.mm(E_torch.view(p1, LibS), out_LR.view(LibS, nr1 * nc1))
        # Smoothing
        if out_avg is None:
            out_avg = out_LR.detach()
            out_HR_avg = out_HR.detach()
        else:
            out_avg = out_avg * exp_weight + out_LR.detach() * (1 - exp_weight)
            out_HR_avg = out_HR_avg * exp_weight + out_HR.detach() * (1 - exp_weight)

        # %%
        out_HR = out_HR.view((1, p1, nr1, nc1))
        total_loss = mse(img_noisy_torch, out_HR)
        total_loss.backward()
        if 1:
            out_LR_np = out_LR.detach().cpu().squeeze().numpy()
            out_avg_np = out_avg.detach().cpu().squeeze().numpy()

            print('Iteration %05d    Loss %f' % (i, total_loss.item()), '\r', end='')
        if PLOT and i % show_every == 0:
            out_LR_np = torch_to_np(out_LR)
            out_avg_np = torch_to_np(out_avg)
            plt.plot(out_LR_np.reshape(LibS, nr1 * nc1))
        loss[i] = total_loss.item()
        i += 1

        return total_loss


    p11 = get_params(OPT_OVER, net1, net_input1)
    optimize(OPTIMIZER1, p11, closure1, LR1, num_iter1)
    if 1:
        out_LR_np = out_LR.detach().cpu().squeeze().numpy()
        out_avg_np = out_avg.detach().cpu().squeeze().numpy()

        print('Iteration %05d' % (i),
              '\r', end='')
    if  save_result is True:
         scipy.io.savemat(r"D:\work\TGRS\data\jasperRidge2_R198\out_avg_np_jasperRidge2_%01d%01d.mat" % (fi+2,1),
                        {'out_avg_np_jasperRidge2_%01d%01d' % (fi+2,1):out_avg_np.transpose(1,2,0)})
         scipy.io.savemat(r"D:\work\TGRS\data\jasperRidge2_R198\out_LR_np_jasperRidge2_%01d%01d.mat" % (fi+2,1),
                        {'out_LR_np_jasperRidge2_%01d%01d' % (fi+2, 1):out_LR_np.transpose(1,2,0)})





"""
用于测试ADSU算法 samson数据

"""


torch.backends.cudnn.enabled = True
torch.backends.cudnn.benchmark =True
dtype = torch.cuda.FloatTensor

PLOT = False
import scipy.io
#%% 加载数据 方块75*75数据(图像数据)
fname2  = "ADSU_Data/real_datasets/samson_1.mat"
mat2 = scipy.io.loadmat(fname2)
img_np_gt = mat2["Y"]
img_np_gt = img_np_gt.transpose(2,0,1)
[p1, nr1, nc1] = img_np_gt.shape

# 加载光谱库
fname4  = "ADSU_Data/real_datasets/EE_samson.mat"
mat4 = scipy.io.loadmat(fname4)
EE = mat4["EE"]

#%% 光谱库中端元数
LibS=EE.shape[1]
save_result = True

import time
from tqdm import tqdm
# rmax = 5;
tol1 = 1

for fi in tqdm(range(tol1)):
    img_noisy_np = img_np_gt
    INPUT = 'noise'  # 'meshgrid'
    pad = 'reflection'
    need_bias = True
    OPT_OVER = 'net'  # 'net,input'

    #
    reg_noise_std = 0.0
    LR1 = 0.001

    OPTIMIZER1 = 'adam'  # 'RMSprop'#'adam' # 'LBFGS'
    show_every = 500
    exp_weight = 0.99

    if fi == 0:
        num_iter1 = 20000
    input_depth =img_noisy_np.shape[0]
    class CAE_AbEst(nn.Module):
        def __init__(self):
            super(CAE_AbEst, self).__init__()
            self.conv1 = nn.Sequential(
                UnmixArch(
                    input_depth, EE.shape[1],
                    num_channels_down=[256],
                    num_channels_up=[256],
                    num_channels_skip=[4],
                    filter_size_up=3, filter_size_down=3, filter_skip_size=1,
                    upsample_mode='bilinear',  # downsample_mode='avg',
                    need1x1_up=True,
                    need_sigmoid=True, need_bias=True, pad=pad, act_fun='LeakyReLU').type(dtype)
            )

        def forward(self, x):
            x = self.conv1(x)
            return x


    net1 = CAE_AbEst()
    net1.cuda()
    print(net1)

    # Compute number of parameters
    s = sum([np.prod(list(p11.size())) for p11 in net1.parameters()]);
    print('Number of params: %d' % s)

    # Loss
    mse = torch.nn.MSELoss().type(dtype)
    img_noisy_torch = np_to_torch(img_noisy_np).type(dtype)
    # if fk==0:
    net_input1 = get_noise(input_depth, INPUT,
                           (img_noisy_np.shape[1], img_noisy_np.shape[2])).type(dtype).detach()
    net_input1 = img_noisy_torch
    E_torch = np_to_torch(EE).type(dtype)
    # %%
    net_input_saved = net_input1.detach().clone()
    noise = net_input1.detach().clone()
    out_avg = None
    out_HR_avg = None
    last_net = None
    RMSE_LR_last = 0
    loss = np.zeros((num_iter1, 1))
    AE = np.zeros((num_iter1, 1))
    i = 0


    def closure1():

        global i, RMSE_LR, RMSE_LR_ave, RMSE_HR, out_LR_np, out_avg_np, out_LR \
            , out_avg, out_HR_np, out_HR_avg, out_HR_avg_np, RMSE_LR_last, last_net \
            , net_input, RMSE_LR_avg, RMSE_HR_avg, RE_HR_avg, RE_HR, Eest, loss, AE \
            , MAE_LR, MAE_LR_avg, MAE_HR, MAE_HR_avg

        if reg_noise_std > 0:
            net_input = net_input_saved + (noise.normal_() * reg_noise_std)

        out_LR = net1(net_input1)
        out_HR = torch.mm(E_torch.view(p1, LibS), out_LR.view(LibS, nr1 * nc1))
        # Smoothing
        if out_avg is None:
            out_avg = out_LR.detach()
            out_HR_avg = out_HR.detach()
        else:
            out_avg = out_avg * exp_weight + out_LR.detach() * (1 - exp_weight)
            out_HR_avg = out_HR_avg * exp_weight + out_HR.detach() * (1 - exp_weight)

        # %%
        out_HR = out_HR.view((1, p1, nr1, nc1))
        total_loss = mse(img_noisy_torch, out_HR)
        total_loss.backward()
        if 1:
            out_LR_np = out_LR.detach().cpu().squeeze().numpy()
            out_avg_np = out_avg.detach().cpu().squeeze().numpy()

            print('Iteration %05d    Loss %f' % (i, total_loss.item()), '\r', end='')
        if PLOT and i % show_every == 0:
            out_LR_np = torch_to_np(out_LR)
            out_avg_np = torch_to_np(out_avg)
            plt.plot(out_LR_np.reshape(LibS, nr1 * nc1))
        loss[i] = total_loss.item()
        i += 1

        return total_loss


    p11 = get_params(OPT_OVER, net1, net_input1)
    optimize(OPTIMIZER1, p11, closure1, LR1, num_iter1)
    if 1:
        out_LR_np = out_LR.detach().cpu().squeeze().numpy()
        out_avg_np = out_avg.detach().cpu().squeeze().numpy()

        print('Iteration %05d' % (i),
              '\r', end='')
    if  save_result is True:
         scipy.io.savemat(r"D:\work\TGRS\data\samson_1\out_avg_np_samson_%01d%01d.mat" % (fi+2,1),
                        {'out_avg_np_samson_%01d%01d' % (fi+2,1):out_avg_np.transpose(1,2,0)})
         scipy.io.savemat(r"D:\work\TGRS\data\samson_1\\out_LR_np_samson_%01d%01d.mat" % (fi+2,1),
                        {'out_LR_np_samson_%01d%01d' % (fi+2, 1):out_LR_np.transpose(1,2,0)})