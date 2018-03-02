import neuron
import pylab as plt
import LFPy
import scipy.io


cell_parameters = {
    'morphology': 'neuron_model/morphologies/cell1.asc',
    # 'delete_sections' : [False]
}

cell = LFPy.Cell(**cell_parameters)

plt.figure()
for idx in range(cell.totnsegs):
    plt.plot([cell.xstart[idx], cell.xend[idx]],
             [cell.ystart[idx], cell.yend[idx]], c='k')
plt.show()

# xstart = cell.xstart
# xend = cell.xend
#
# ystart = cell.ystart
# yend = cell.yend
#
# zstart = cell.zstart
# zend = cell.zend
#
# N = cell.totnsegs


savedict = {
    'xstart': cell.xstart,
    'xend': cell.xend,
    'ystart': cell.ystart,
    'yend': cell.yend,
    'zstart': cell.zstart,
    'zend': cell.zend,
    'N': cell.totnsegs
}

# scipy.io.savemat('cell_morphology.mat', savedict)
