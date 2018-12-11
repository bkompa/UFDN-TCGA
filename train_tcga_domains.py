import yaml
import os
import sys
import shutil
import numpy as np
import torch

from itertools import cycle

from torch.backends import cudnn
import torch.optim as optim
import torch.nn as nn
from torch.autograd import Variable, grad

from src.data import LoadDataset
from src.ufdn import LoadModel

#We'll need to modify this for 3+ domains
from src.util import vae_loss, calc_gradient_penalty, interpolate_vae_3d

from tensorboardX import SummaryWriter 


# Experiment Setting
cudnn.benchmark = True
config_path = sys.argv[1]
conf = yaml.load(open(config_path,'r'))
exp_name = conf['exp_setting']['exp_name']
#img_size is only used in conv nets
#originally it was 64
img_size = conf['exp_setting']['img_size']
#20,501 img_depth
img_depth = conf['exp_setting']['img_depth']

#dictionary with number of steps and batch size
trainer_conf = conf['trainer']

if trainer_conf['save_checkpoint']:
    model_path = conf['exp_setting']['checkpoint_dir'] + exp_name+'/'
    if not os.path.exists(model_path):
        os.makedirs(model_path)
    model_path = model_path+'{}'

if trainer_conf['save_log'] or trainer_conf['save_fig']:
    if os.path.exists(conf['exp_setting']['log_dir']+exp_name):
        shutil.rmtree(conf['exp_setting']['log_dir']+exp_name)
    writer = SummaryWriter(conf['exp_setting']['log_dir']+exp_name)


# Fix seed
np.random.seed(conf['exp_setting']['seed'])
_ = torch.manual_seed(conf['exp_setting']['seed'])

# Load dataset

#TODO: add more than 3 domains 
#text of domain name used for h5 file
#conf['exp_setting']['domains']
#domains: ['brca', 'dlbc',...,'acc']
domains = conf['exp_setting']['domains']
number_of_domains = len(domains)


data_root = conf['exp_setting']['data_root']
batch_size = conf['trainer']['batch_size']

#TODO: add more domains 
#list of loaders? 
#loaders = [Load for cancer in domains]

train_loaders = [LoadDataset('tcga', data_root, batch_size, 'train', style=domains[0])]
train_loaders += [cycle(LoadDataset('tcga',data_root,batch_size,'train',style=domain)) for domain in domains[1:]] 


test_loaders = [LoadDataset('tcga',data_root,batch_size,'test',style=domain) for domain in domains] 


# Load Model
enc_dim = conf['model']['vae']['encoder'][-1][1] #latent space dimension #100
code_dim = conf['model']['vae']['code_dim'] #number of domains #currently 3 
vae_learning_rate = conf['model']['vae']['lr'] #learning rate #10e-4
vae_betas = tuple(conf['model']['vae']['betas']) #used for adam optimizer
df_learning_rate = conf['model']['D_feat']['lr'] #feature discriminator LR 
df_betas = tuple(conf['model']['D_feat']['betas'])
dp_learning_rate = conf['model']['D_pix']['lr']
dp_betas = tuple(conf['model']['D_pix']['betas'])

vae = LoadModel('vae',conf['model']['vae'],img_size,img_depth) #img_size=64, img_depth = 3. Last arg is input dimension
d_feat = LoadModel('nn',conf['model']['D_feat'],img_size,enc_dim) #img_size=64, enc_dim=1024
d_pix = LoadModel('nn',conf['model']['D_pix'],img_size,img_depth) #img_size=64, img_depth = 3

reconstruct_loss = torch.nn.MSELoss()
clf_loss = nn.BCEWithLogitsLoss()


# Use cuda
vae = vae.cuda()
d_feat = d_feat.cuda()
d_pix = d_pix.cuda()

reconstruct_loss = reconstruct_loss.cuda()
clf_loss = clf_loss.cuda()


# Optmizer
opt_vae = optim.Adam(list(vae.parameters()), lr=vae_learning_rate, betas=vae_betas)
opt_df = optim.Adam(list(d_feat.parameters()), lr=df_learning_rate, betas=df_betas)
opt_dp = optim.Adam(list(d_pix.parameters()), lr=dp_learning_rate, betas=dp_betas)

# Training

vae.train()
d_feat.train()
d_pix.train()

    
# Domain code setting
domain_code = np.identity(code_dim)

domain_code = torch.FloatTensor(domain_code)



# with more domains, this 





# Loss weight setting
loss_lambda = {}
for k in trainer_conf['lambda'].keys():
    init = trainer_conf['lambda'][k]['init']
    final = trainer_conf['lambda'][k]['final']
    step = trainer_conf['lambda'][k]['step']
    loss_lambda[k] = {}
    loss_lambda[k]['cur'] = init
    loss_lambda[k]['inc'] = (final-init)/step
    loss_lambda[k]['final'] = final



# Training 
global_step = 0

while global_step < trainer_conf['total_step']: 
    if global_step%2 == 0: 
      perm = np.random.permutation(code_dim)

      forword_code = np.identity(code_dim)[perm, :]
      backword_code = forword_code.T 

      forword_code = torch.FloatTensor(forword_code)
      backword_code = torch.FloatTensor(backword_code)

    for train_images in zip(*train_loaders):
        
        # data augmentation
        input_img = torch.cat([image.type(torch.FloatTensor) for image in train_images],dim=0)
        input_img =  Variable(input_img.cuda(),requires_grad=False)


        code = Variable(torch.FloatTensor(domain_code).cuda(),requires_grad=False)
        invert_code = 1-code

        if global_step%2 == 0:
            trans_code = Variable(torch.FloatTensor(forword_code).cuda(),requires_grad=False)
        else:
            trans_code = Variable(torch.FloatTensor(backword_code).cuda(),requires_grad=False) 
      
  
        # Train Feature Discriminator
        opt_df.zero_grad()
        
        enc_x = vae(input_img,return_enc=True).detach()
        code_pred = d_feat(enc_x)

        df_loss = clf_loss(code_pred,code)
        df_loss.backward()
     
        opt_df.step()
        
        # Train Pixel Discriminator
        opt_dp.zero_grad()
        
        pix_real_pred,pix_real_code_pred = d_pix(input_img)
        
        fake_img = vae(input_img,insert_attrs=trans_code)[0].detach()
        pix_fake_pred, _  = d_pix(fake_img)
        
        pix_real_pred = pix_real_pred.mean()
        pix_fake_pred = pix_fake_pred.mean()

        gp = loss_lambda['gp']['cur']*calc_gradient_penalty(d_pix,input_img.data,fake_img.data)
        pix_code_loss = clf_loss(pix_real_code_pred,code)
        
        d_pix_loss = pix_code_loss + pix_fake_pred - pix_real_pred + gp
        d_pix_loss.backward()
        
        opt_dp.step()
        

        # Train VAE
        opt_vae.zero_grad()
        
        ### Reconstruction Phase
        recon_batch, mu, logvar = vae(input_img,insert_attrs = code)
        mse,kl = vae_loss(recon_batch, input_img, mu, logvar, reconstruct_loss)  #.view(batch_size,-1)
        recon_loss = (loss_lambda['pix_recon']['cur']*mse+loss_lambda['kl']['cur']*kl)
        recon_loss.backward()

        
        ### Feature space adversarial Phase       
        enc_x = vae(input_img,return_enc=True)
        domain_pred = d_feat(enc_x)
        adv_code_loss = clf_loss(domain_pred,invert_code)
        
        feature_loss = loss_lambda['feat_domain']['cur']*adv_code_loss
        feature_loss.backward()
        
        ### Pixel space adversarial Phase
        enc_x = vae(input_img,return_enc=True).detach()
        
        fake_img = vae.decode(enc_x,trans_code)
        recon_enc_x = vae(fake_img,return_enc=True)
        adv_pix_loss, pix_code_pred = d_pix(fake_img)
        adv_pix_loss = adv_pix_loss.mean()
        pix_clf_loss = clf_loss(pix_code_pred,trans_code)
        
        
        pixel_loss =  - loss_lambda['pix_adv']['cur']*adv_pix_loss + loss_lambda['pix_clf']['cur']*pix_clf_loss
        pixel_loss.backward()
        
        opt_vae.step()
        
        
        # End of step      
        print('Step',global_step,end='\r',flush=True)     
        global_step += 1
        
        # Records
        if trainer_conf['save_log'] and (global_step % trainer_conf['verbose_step'] ==0):
            writer.add_scalar('MSE', mse.data[0], global_step)
            writer.add_scalar('KL',  kl.data[0], global_step)
            writer.add_scalar('gp', gp.data[0], global_step)
            writer.add_scalars('Pixel_Distance',{'real':pix_real_pred.data[0],
                                               'fake':pix_fake_pred.data[0]}, global_step)
            writer.add_scalars('Code_loss',{'feature':df_loss.data[0],
                                            'pixel':pix_code_loss.data[0],
                                           'adv_feature':feature_loss.data[0],
                                           'adv_pixel':pix_clf_loss.data[0]}, global_step)

            
        # update lambda
        for k in loss_lambda.keys():
            if loss_lambda[k]['inc']*loss_lambda[k]['cur'] < loss_lambda[k]['inc']*loss_lambda[k]['final']:
                loss_lambda[k]['cur'] += loss_lambda[k]['inc']
                


        if global_step%trainer_conf['checkpoint_step']==0 and trainer_conf['save_checkpoint'] and not trainer_conf['save_best_only']:
            torch.save(vae,model_path.format(global_step)+'.vae')

