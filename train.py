import os
import os.path as osp
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from datetime import datetime
from tqdm import tqdm
import argparse
from pathlib import Path
import numpy as np
import zarr
from pathlib import Path

from utils.cov import coverage_score as cov

from utils.utils import seed_everything, Wandb_logger, _optimizer, coverage_score
from model.pipeline import Pipeline
from model.graph_gmvae import DeepMetaBinModel
import shutil

def main():
    ######Hyperparameters######
    parser = argparse.ArgumentParser()
    # model
    parser.add_argument("--GPU", default=0, type=int, help="GPU index")

    parser.add_argument("--KNN", "-k", default=3, type=int, help="K-nearest neighbor")

    parser.add_argument("--learning_rate", "-lr", default=1e-4, type=float, help="Learning rate")
    parser.add_argument("--weight_decay", default=5e-4, type=float, help="Weight decay for optimizer")
    parser.add_argument("--w_cat", default=1, type=float, help="Weight for Categorical loss")
    parser.add_argument("--w_gauss", default=1, type=float, help="Weight for Gaussian loss")
    parser.add_argument("--w_rec", default=1, type=float, help="Weight for Reconstruction loss")
    parser.add_argument("--w_cl", default=0.5, type=float, help="Weight for Reconstruction loss")
    parser.add_argument("--input_size", default=104, type=int, help="Input feature size")
    parser.add_argument("--gaussian_size", default=32, type=int, help="Embed size")
    parser.add_argument("--sigma", default=1.0, type=float, help="The sigma for Gassian kernal")

    #data
    parser.add_argument("--seed", type=int, default=2024, help="Seed")
    parser.add_argument("--wandb", type=str, default='disabled', choices=['online', 'offline', 'disabled', 'dryrun'])
    parser.add_argument("--zarr_dataset_path", '-data', type=str, default='', help="Dataset zarr path")
    parser.add_argument("--contignames_path", type=str, default='./sample_data/contignames.npz', help="Contigname path")
    parser.add_argument("--contig_path", type=str, default='./sample_data/contigs.fasta', help="Contig fasta path")
    parser.add_argument("--exp_name", "-exp", type=str, default='time', help="Name for this experiment")
    parser.add_argument("--batch_size", "-b", type=int, default=420, help="Batch size for NN")
    parser.add_argument("--num_workers", type=int, default=50, help="Number of workers")
    parser.add_argument("--output", type=str, default="./deepmetabin_out", help="Output for deepmetabin")
    parser.add_argument("--num_epoch", "-e", type=int, default=500, help="Epoch for NN")
    parser.add_argument("--multisample", type=bool, default=False, help="Multi-sample or single-sample")
    args = parser.parse_args()

    ######Initialization######
    if args.exp_name == 'time':
        args.exp_name = datetime.now().strftime("%Y%m%d%H%M%S")
    working_dir = os.path.join('./exp', args.exp_name)
    Path(working_dir).mkdir(parents=True, exist_ok=True)
    logging = Wandb_logger(cfg=vars(args), working_dir=working_dir)
    if torch.cuda.is_available():
        gpu_index = args.GPU
        device = f"cuda:{gpu_index}"
    else:
        device = "cpu"
    seed_everything(args.seed)

    ######Pipeline######

    os.makedirs(args.output, exist_ok=True)
    pip = Pipeline(zarr_dataset_path=args.zarr_dataset_path,
                   k=args.KNN,
                   sigma=args.sigma,
                   multisample=args.multisample,
                   must_link_path=osp.join(args.output, 'must_link.csv')
                   )
    dataloader = DataLoader(
                    dataset=pip,
                    batch_size=args.batch_size,
                    num_workers=args.num_workers,
                    pin_memory=False,
                    shuffle=True,
                    )
    val_loader = DataLoader(
                    dataset=pip,
                    batch_size=len(pip),
                    num_workers=args.num_workers,
                    pin_memory=False,
                    shuffle=False,
                    )

    model = DeepMetaBinModel(input_size=args.input_size,
                             gaussian_size=args.gaussian_size,
                             w_cat=args.w_cat,
                             w_gauss=args.w_gauss,
                             w_rec=args.w_rec,
                             w_cl=args.w_cl,
                             zarr_dataset_path=args.zarr_dataset_path,
                             contignames_path=args.contignames_path,
                             log_path=args.output,
                            #  use_gmm=True
                             k=args.KNN,
                             result_path=osp.join(args.output, 'results'),
                             contig_path=args.contig_path
                             )
    scheduler, optimizer = _optimizer(model=model, 
                        lr=args.learning_rate, 
                        weight_decay=args.weight_decay, 
                        epoch=args.num_epoch)
     ######Training the model######
    logging.info("Start Training...")

    for epoch in range(args.num_epoch):
        logging.info(f"Epoch ({epoch}/{args.num_epoch})")
        model.train()
        for i, batch in enumerate(tqdm(dataloader, ncols=80, desc='Training')):
            optimizer.zero_grad()
            loss = model.training_step(batch, i)['loss']
            logging.logging_with_step('loss', loss, epoch * len(dataloader) + i)
            loss.backward()
            optimizer.step()
        logging.info(f'loss: {loss}')
        scheduler.step()

        
    model.eval()
    model.validation_step()
    logging.info("Wrote contigs into bins")
    logging.info('Finish training!')

    # Training include early stopping(deactivated)
    # patience_counter = 0
    # best_coverage = 0
    # patience = 5
    # epoch_c = 0
    # gmmcsv_path = os.path.join(args.output, 'results/gmm.csv')
    # gmmcsv_best_path = os.path.join(args.output, 'results/gmm_best.csv')
    # metrics = []
    # for epoch in range(args.num_epoch):
    #     logging.info(f"Epoch ({epoch}/{args.num_epoch})")
    #     model.train()
    #     for i, batch in enumerate(tqdm(dataloader, ncols=80, desc='Training')):
    #         optimizer.zero_grad()
    #         loss = model.training_step(batch, i)['loss']
    #         logging.logging_with_step('loss', loss, epoch * len(dataloader) + i)
    #         loss.backward()
    #         optimizer.step()
    #     logging.info(f'loss: {loss}')
    #     scheduler.step()

    #     model.eval()
    #     if epoch_c >= 100 and epoch_c % 5 == 0:
    #         with torch.no_grad():
    #             for batch in val_loader:
    #                 model.validation_step(batch)
    #                 metrics.append(cov(gmmcsv_path))
            
    #         coverage = coverage_score(gmmcsv_path)

    #         if coverage >= best_coverage:
    #             best_coverage = coverage
    #             patience_counter = 0
    #             torch.save(model.state_dict(), 'best_model.pth')
    #             shutil.copy(gmmcsv_path, gmmcsv_best_path)
    #         else:
    #             patience_counter += 1
    #             if patience_counter >= patience:
    #                 model.load_state_dict(torch.load('best_model.pth'))
                    
    #                 model.eval()
    #                 with torch.no_grad():
    #                     for batch in val_loader:
    #                         model.validation_step(batch)
                    
    #                 logging.info('Early stopping triggered, best model restored')
    #                 shutil.copy(gmmcsv_best_path, gmmcsv_path)
    #                 model.rec_best_gmm()

    #                 metrics.append(cov(gmmcsv_path))
    #                 # Write metrics into a CSV file
    #                 metrics_path = os.path.join(args.output, 'metrics.csv')
    #                 metrics_array = np.array(metrics)
    #                 np.savetxt(metrics_path, metrics_array, delimiter=',', fmt='%1.8f')

    #                 break
    #     epoch_c += 1
        # logging.info("Wrote contigs into bins")
        # logging.info('Finish training!')
    

    



if __name__ == '__main__':
    main()