import os
from pathlib import Path
from pfs.drp.stella.referenceLine import ReferenceLineSet

ref_linelist_dir = os.path.join(os.environ['OBS_PFS_DIR'], 'pfs', 'lineLists')
fileList = ['Hg', 'Ar']

lines = []
for filename in fileList:
    for r in ReferenceLineSet.fromLineList(Path(ref_linelist_dir,
                                                filename).with_suffix('.txt')):
        lines.append(r)

print(type(lines[0]))
linesSorted = sorted(lines, key=lambda w: w.wavelength)

referencelineSet = ReferenceLineSet.empty()

for r in linesSorted:
    referencelineSet.append(description=r.description,
                            wavelength=r.wavelength,
                            intensity=r.intensity,
                            status=r.status)

referencelineSet.writeLineList('out.txt')
