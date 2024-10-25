import m2aia as m2
import numpy as np
import tensorflow as tf
import SimpleITK as sitk

from msiPL.Computational_Model import VAE_BN
from msiPL.LearnPeaks import *

from pathlib import Path
import argparse

# create a BatchSequence for keras using pym2aia's BatchGenerator
class BatchSequence(tf.keras.utils.Sequence):
    def __init__(self, dataset: m2.Dataset.BaseDataSet , batch_size: int, shuffle: bool=True):
        super().__init__()
        self.gen = m2.BatchGenerator(dataset, batch_size, shuffle)
    
    def __len__(self):
        return self.gen.__len__()

    def on_epoch_end(self):
        self.gen.on_epoch_end()
    
    def __getitem__(self, index):
        return self.gen.__getitem__(index)[0]

def running_variance_update(existingAggregate, newValue):
    (count, mean, deltaM2) = existingAggregate
    count += 1
    delta = newValue - mean
    mean += delta / count
    delta2 = newValue - mean
    deltaM2 += delta * delta2
    return (count, mean, deltaM2)

def running_variance_finalize(existingAggregate):
    (count, mean, deltaM2) = existingAggregate
    (mean, variance, sampleVariance) = (mean, deltaM2 / count, deltaM2 / (count - 1))
    return (mean, variance, sampleVariance)


# def add_qt_widget(func):
#     def wrapper(*args, **kwargs):
#         print(args, kwargs)
#         return func(*args, **kwargs)
#     return wrapper


parser = argparse.ArgumentParser(prog="msiPL", 
                                 description="Running msiPL!",
                                 epilog="Original work by Abdelmoula, Walid M. et al. Peak Learning of Mass Spectrometry Imaging Data Using Artificial Neural Networks. Nature Communications 12, no. 1 (December 2021): 5544. https://doi.org/10.1038/s41467-021-25744-8.")

# parser.add_argument = add_qt_widget(parser.add_argument)

parser.add_argument("-i", "--imzml", required=True)
parser.add_argument("--mask", required=True)
parser.add_argument("-o", "--csv", required=True)
parser.add_argument("--model", default="/tmp/msiPL.h5", required=False)
parser.add_argument("-l", "--latent_dim", default=5, required=False, type=int)
parser.add_argument("-d", "--interim_dim", default=256, required=False, type=int)
parser.add_argument("-b", "--beta", default=1.3, required=False, type=float)
parser.add_argument("-e", "--epochs", default=2, required=False, type=int)
parser.add_argument("-s", "--batch_size", default=128, required=False, type=int)

parser.add_argument("--normalization", default="None")
parser.add_argument("--intensity_transform", default="None" )
parser.add_argument("--range_pooling", default="Mean")
parser.add_argument("--smoothing", default="None")
parser.add_argument("--smoothing_value", default=2, type=np.int32)
parser.add_argument("--baseline_correction", default="None")
parser.add_argument("--baseline_correction_value", default=50, type=np.int32)
parser.add_argument("--tolerance", default=75, type=np.float32)

args = parser.parse_args()

latent_dim = args.latent_dim
interim_dim = args.interim_dim
epochs = args.epochs
batch_size = args.batch_size
file_path : str = args.imzml
beta = args.beta 

model_path = args.model

itkMask = sitk.ReadImage(args.mask)
mask = sitk.GetArrayFromImage(itkMask)

# load data
I = m2.ImzMLReader(file_path)
I.SetSmoothing(args.smoothing, args.smoothing_value)
I.SetBaselineCorrection(args.baseline_correction, args.baseline_correction_value)
I.SetNormalization(args.normalization)
I.SetPooling(args.range_pooling)
I.SetIntensityTransformation(args.intensity_transform)
I.SetTolerance(args.tolerance)
I.Load()


# check for existing model parameters
file_name = Path(file_path).name

print(mask.shape, (mask == 1).shape)
print(np.unique(mask))

for label in np.unique(mask):
    if label == 0: 
        continue
    
    # initialize model
    vae = VAE_BN(I.GetXAxisDepth(), interim_dim, latent_dim)
    myModel, encoder = vae.get_architecture()
    
    print(f"Training model for label {label}...")
    print(f"Model will be saved at {model_path}")
    print("Samples in mask:", np.sum(mask == label))

    dataset = m2.Dataset.SpectrumDataset([I],sampling_masks=[mask == label])
    gen = BatchSequence(dataset, batch_size=batch_size, shuffle=True)
    history = myModel.fit(gen, epochs=epochs)
    myModel.save_weights(model_path)

    # init running variance
    count = 0
    mean = np.zeros_like(I.GetXAxis())
    deltaM2  = np.zeros_like(I.GetXAxis())
    existingAggregate = (count, mean, deltaM2)

    # update running variance
    for i in range(I.GetNumberOfSpectra()):
        xs, ys = I.GetSpectrum(i)
        existingAggregate = running_variance_update(existingAggregate, ys)

    _, var, _ = running_variance_finalize(existingAggregate)

    # evaluate encoder weights and find peaks
    W_enc = encoder.get_weights()
    _, _, _, learned_peaks = LearnPeaks(I.GetXAxis(), W_enc, np.sqrt(var), latent_dim, beta, I.GetMeanSpectrum())

    # write as csv
    xs = I.GetXAxis()
    ys = I.GetMeanSpectrum()
    name = args.csv.split('.')[0] + "_" + str(label) + ".csv"

    with open(name, 'w') as f:
        f.write('mz,intensity\n')
        for id in learned_peaks:
            f.write(f'{xs[id]},{ys[id]}\n')