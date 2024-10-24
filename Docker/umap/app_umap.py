import argparse
import numpy as np
parser = argparse.ArgumentParser(prog="UMAP", 
                                 description="Running simCLR!",
                                 epilog="https://github.com/lmcinnes/umap/blob/master/umap/umap_.py#L1403")

# parser.add_argument = add_qt_widget(parser.add_argument)

parser.add_argument("-i", "--imzml", required=True)

parser.add_argument("--normalization", default="None")
parser.add_argument("--intensity_transform", default="None" )
parser.add_argument("--range_pooling", default="Mean")
parser.add_argument("--smoothing", default="None")
parser.add_argument("--baseline_correction", default="None")
parser.add_argument("--smoothing_value", default=2, type=np.int32)
parser.add_argument("--baseline_correction_value", default=50, type=np.int32)
parser.add_argument("--tolerance", default=75, type=np.float32)

parser.add_argument("-c", "--centroids", required=True)
parser.add_argument("-o", "--out", required=True)

parser.add_argument("-e","--n_epochs", default=None, type=int)
parser.add_argument("-nn","--n_neighbors", default=15, type=np.int32)
parser.add_argument("-nc","--n_components", default=2, type=np.int32)
parser.add_argument("-j","--n_jobs", default=-1, type=np.int32)
parser.add_argument("-m","--metric", default="euclidean", type=str)
parser.add_argument("-mkwds","--metric_kwds", default=None, type=str)
parser.add_argument("-om","--output_metric", default="euclidean", type=str)
parser.add_argument("-omkwds","--output_metric_kwds", default=None, type=str)
parser.add_argument("-lr","--learning_rate", default=1.0, type=np.float32)
parser.add_argument("-in","--init", default="spectral", type=str)
parser.add_argument("-md","--min_dist", default=0.1, type=np.float32)
parser.add_argument("-s","--spread", default=1.0, type=np.float32)
parser.add_argument("-lm","--low_memory", default=True, action="store_true")
parser.add_argument("--set_op_mix_ratio", default=1.0, type=np.float32)
parser.add_argument("--local_connectivity", default=1.0, type=np.float32)
parser.add_argument("--repulsion_strength", default=1.0, type=np.float32)
parser.add_argument("--negative_sample_rate", default=5, type=np.float32)
parser.add_argument("--transform_queue_size", default=4.0, type=np.float32)
parser.add_argument("--a", default=None, type=np.float32)
parser.add_argument("--b", default=None, type=np.float32)
parser.add_argument("-r","--random_state", default=None, type=np.int32)
parser.add_argument("--angular_rp_forest", default=False, action="store_true")
parser.add_argument("--target_n_neighbors", default=-1, type=np.int32)
parser.add_argument("--target_metric", default="categorical", type=str)
parser.add_argument("--target_metric_kwds", default=None, type=str)
parser.add_argument("--target_weight", default=0.5, type=np.float32)
parser.add_argument("--transform_seed", default=42, type=np.int32)
parser.add_argument("--transform_mode", default="embedding", type=str)
parser.add_argument("--force_approximation_algorithm", default=False, action="store_true")
parser.add_argument("--verbose", default=False, action="store_true")
parser.add_argument("--tqdm_kwds", default=None, type=str)
parser.add_argument("--unique", default=False, action="store_true")
parser.add_argument("-d","--densmap", default=False, action="store_true")
parser.add_argument("-dl","--dens_lambda", default=2.0, type=np.float32)
parser.add_argument("-df","--dens_frac", default=0.3, type=np.float32)
parser.add_argument("-dvs","--dens_var_shift", default=0.1, type=np.float32)
parser.add_argument("-do","--output_dens", default=False, action="store_true")
parser.add_argument("--disconnection_distance", default=None, type=np.float32)

args = parser.parse_args()

print(args, type(args.n_epochs))

import m2aia as m2
import umap

model = umap.UMAP(
        n_neighbors=args.n_neighbors,
        n_components=args.n_components,
        metric=args.metric,
        metric_kwds=args.metric_kwds,
        output_metric=args.output_metric,
        output_metric_kwds=args.output_metric_kwds,
        n_epochs=args.n_epochs,
        learning_rate=args.learning_rate,
        init=args.init,
        min_dist=args.min_dist,
        spread=args.spread,
        low_memory=args.low_memory,
        n_jobs=args.n_jobs,
        set_op_mix_ratio=args.set_op_mix_ratio,
        local_connectivity=args.local_connectivity,
        repulsion_strength=args.repulsion_strength,
        negative_sample_rate=args.negative_sample_rate,
        transform_queue_size=args.transform_queue_size,
        a=args.a,
        b=args.b,
        random_state=args.random_state,
        angular_rp_forest=args.angular_rp_forest,
        target_n_neighbors=args.target_n_neighbors,
        target_metric=args.target_metric,
        target_metric_kwds=args.target_metric_kwds,
        target_weight=args.target_weight,
        transform_seed=args.transform_seed,
        transform_mode=args.transform_mode,
        force_approximation_algorithm=args.force_approximation_algorithm,
        verbose=args.verbose,
        tqdm_kwds=args.tqdm_kwds,
        unique=args.unique,
        densmap=args.densmap,
        dens_lambda=args.dens_lambda,
        dens_frac=args.dens_frac,
        dens_var_shift=args.dens_var_shift,
        output_dens=args.output_dens,
        disconnection_distance=args.disconnection_distance)


print("Load: ", args.imzml)
I = m2.ImzMLReader(args.imzml)
I.SetSmoothing(args.smoothing, args.smoothing_value)
I.SetBaselineCorrection(args.baseline_correction, args.baseline_correction_value)
I.SetNormalization(args.normalization)
I.SetPooling(args.range_pooling)
I.SetIntensityTransformation(args.intensity_transform)
I.Load()


# Load centroids
C = np.genfromtxt(args.centroids, dtype=np.float32, delimiter=',', skip_header=1)

# Get the mask
M = I.GetMaskArray()
N = np.sum(M)

D = np.zeros((C.shape[0], N))
for i,c in enumerate(C[:,0]):
    D[i] = I.GetArray(c, args.tolerance)[M==1]

# apply umap to data
Y = model.fit_transform(D.transpose())
s = M.shape
O = np.zeros((s[0], s[1], s[2], args.n_components))

# reshape data back to image sapce
for i in range(args.n_components):
    Y[...,i] = Y[...,i]-np.min(Y[...,i])
    Y[...,i] = Y[...,i]/np.max(Y[...,i])*255.0

    O[...,i][M==1] = Y[...,i]

# write the image as a vector image
import SimpleITK as sitk
OI = sitk.GetImageFromArray(O, True)
sitk.WriteImage(OI, args.out)
