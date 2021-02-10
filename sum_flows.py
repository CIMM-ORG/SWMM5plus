import os
import numpy as np

fpath = '/home/griano/Documents/Github/SWMMengine/swmm_files/brazos'

sum_check = 0
norm_infty = 0
for i, p in enumerate(os.listdir(fpath)):
    if ".dat" in p:
        file_path = os.path.join(fpath, p)
        data = np.loadtxt(file_path)
        sum_check += sum(data[:,1])
        mmax = max(data[:,1])
        if norm_infty < mmax:
            norm_infty = mmax
        print(i, p, norm_infty, sum_check)

print("The sum is", sum_check)
