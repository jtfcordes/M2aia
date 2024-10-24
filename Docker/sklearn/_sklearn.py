import m2aia as m2
import numpy as np
import SimpleITK as sitk
from pathlib import Path
from sklearn.cluster import BisectingKMeans, KMeans
from skimage.filters import gaussian
import argparse
import ast

parser = argparse.ArgumentParser(prog="sklearn", 
                                 description="Running sklearn (clustering)!",
                                 epilog="First prototype to run from sklearn.cluster BisectingKMeans and KMeans")

# parser.add_argument = add_qt_widget(parser.add_argument)

parser.add_argument("-i", "--imzml", required=True)
parser.add_argument("-c", "--centroids", required=True)
parser.add_argument("-o", "--out", required=True)
# parser.add_argument("--out_cluster_centers")
# parser.add_argument("--in_cluster_centers")
parser.add_argument("--masking", action='store_true')
parser.add_argument("--gaussian", action='store_true')
parser.add_argument("--sigma", type=lambda x: ast.literal_eval(x), default="[0.008,0.008]")

# clustering parameters
parser.add_argument("--clustering", action='store_true')
parser.add_argument("-n", "--num_clusters")
parser.add_argument("--clustering-method")

parser.add_argument("--normalize", action='store_true')
parser.add_argument("--write-ion-images", action='store_true')
# co-localization
# parser.add_argument("-n", "--num_clusters")
# parser.add_argument("--clustering-method")

parser.add_argument("--normalization", default="None")
parser.add_argument("--intensity_transform", default="None" )
parser.add_argument("--range_pooling", default="Maximum")
parser.add_argument("--smoothing", default="None")
parser.add_argument("--smoothing_value", default=2, type=np.int32)
parser.add_argument("--baseline_correction", default="None")
parser.add_argument("--baseline_correction_value", default=50, type=np.int32)
parser.add_argument("--tolerance", default=0.03, type=np.float32)

args = parser.parse_args()

print("[imzML path]", args.imzml)
print("normalization ->".rjust(30), args.normalization)
print("tolerance ->".rjust(30), args.tolerance)
print("intensity_transform ->".rjust(30), args.intensity_transform)
print("range_pooling ->".rjust(30), args.range_pooling)
print("smoothing ->".rjust(30), args.smoothing)
print("smoothing_value ->".rjust(30), args.smoothing_value)
print("baseline_correction ->".rjust(30), args.baseline_correction)
print("baseline_correction_value ->".rjust(30), args.baseline_correction_value)




print("[centroids path]", args.centroids)
if args.centroids:
    C = np.genfromtxt(args.centroids, dtype=np.float32, delimiter=',', skip_header=1)
    print("number ->".rjust(30), C.shape[0])
    if C.shape[0] < 25:
        print("centroids ->".rjust(30), C[:,0])

print("* [gaussian smoothing]", args.gaussian)
if args.smoothing:
    print("sigma ->".rjust(30), args.sigma)

print("* [clustering]", args.clustering)
if args.clustering:
    print("num_clusters ->".rjust(30), args.num_clusters)

print("* [normalize]", args.normalize)
print("* [write ion-images]", args.write_ion_images)

I = m2.ImzMLReader(args.imzml)
I.SetSmoothing(args.smoothing, args.smoothing_value)
I.SetBaselineCorrection(args.baseline_correction, args.baseline_correction_value)
I.SetNormalization(args.normalization)
I.SetPooling(args.range_pooling)
I.SetIntensityTransformation(args.intensity_transform)
I.Load()

M = I.GetMaskArray()
print("* [external masking]", args.masking)
if args.masking:
    maskPath = args.imzml.split(".imzML")[0] + ".mask.nrrd"
    if Path(maskPath).exists():
        print("mask path ->".rjust(30), maskPath)
        M = sitk.GetArrayFromImage(sitk.ReadImage(maskPath))
N = np.sum(M)
print("masked pixels ->".rjust(30), N)

print("------------------------")
exit(0)

D = np.zeros((C.shape[0] , N))
for i,c in enumerate(C[:,0]):
    A = I.GetImage(c, args.tolerance)
    if args.gaussian:
        A[...,0] = sitk.SmoothingRecursiveGaussian(A[...,0], [0.008,0.008])
    if args.normalize:
        A = sitk.Normalize(A)
    if args.write_ion_images:
        sitk.WriteImage(I.GetImage(c,args.tolerance), Path(args.out).parent/f"ion_{i}_{c}.nrrd")
    D[i] = sitk.GetArrayFromImage(A)[M>=1]
    


# print("row_indices", row_indices.shape)
# print("col_indices", col_indices.shape)
print("D.shape", D.shape)
print("M.shape", M.shape)

# D[-2] = row_indices[M[0]>=1]
# D[-1] = col_indices[M[0]>=1]

algo = BisectingKMeans(n_clusters=int(args.num_clusters), n_init=2, 
bisecting_strategy='largest_cluster')
# algo.cluster_centers_ = [0,0,0]
# print(algo.cluster_centers_)
Y = algo.fit_predict(D.T)

M[M==1] = Y+1

M = sitk.GetImageFromArray(M)
M.CopyInformation(I.GetMaskImage())
sitk.WriteImage(M, args.out)
