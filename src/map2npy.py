import sys
import numpy as np
from chimerax.core.commands import run

target_filename = sys.argv[1]

target_map = run(session, 'open %s' % target_filename)[0]

target_vol = target_map.data.full_matrix()
target_vol_npy_filename = target_filename.rsplit('.', 1)[0] + '.npy'
np.save(target_vol_npy_filename, target_vol)

print("Saved: \n%s\n%s\n%s", target_vol_npy_filename)

