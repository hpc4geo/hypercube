
import numpy as np
import hcell as hyper
import branin as branin
from scipy.interpolate import RBFInterpolator


def point_in_cell(points, tree):
  leaves = hyper.HCell.collect_leaves(tree)
  nl = len(leaves)
  bb0 = np.zeros((nl, 2))
  bb1 = np.zeros((nl, 2))
  for k in range(nl):
    bbmi,bbma = leaves[k].get_bounding_box()
    bb0[k, :] = bbmi[:]
    bb1[k, :] = bbma[:]

  containingcells = list()
  npoints = points.shape[0]
  for p in range(npoints):
    cc = points[p, :]
    found = 0
    found_idx = -1
    for k in range(nl):
      if cc[0] > bb0[k, 0] and cc[1] > bb0[k, 1]:
        if cc[0] < bb1[k, 0] and cc[1] < bb1[k, 1]:
          found += 1
          found_idx = k
    if found > 0:
      if found == 1:
        containingcells.append( leaves[found_idx] )
      else:
        raise RuntimeError('Multiple cells claimed point')


  return containingcells

def t1():

  np.random.seed(0)
  ntest = 100000
  xy_test = np.random.random((ntest, 2))
  xy_test[:, 0] *= 15.0
  xy_test[:, 0] += -5.0
  xy_test[:, 1] *= 15.0
  F_test = branin.evaluate(xy_test[:, 0], xy_test[:, 1])


  root = hyper.HCell([-5.0, 0.0], width=[15.0, 15.0])

  root.refine(nsub=[32, 32])
  cells = root.get_leaves()

  ncells = len(cells)
  #xy = np.zeros((ncells, 2))
  xy = hyper.HCell.get_centroids(cells)

  F = branin.evaluate(xy[:, 0], xy[:, 1])

  rb = RBFInterpolator(xy, F, kernel="quintic")

  #rb(np.array([[1.0, 1.0]]))

  F_est = rb(xy_test)
  error = F_est - F_test
  error = np.absolute(error)
  print(np.max(error))

  idx = np.argsort(error)
  print('err', error[idx[0]], error[idx[-1]])
  print('coor', idx[-1])


  target_xy = xy_test[ idx[-1], : ]
  #target_cell = point_in_cell(np.array([[1,1]]), cells)
  target_cell = point_in_cell(np.array([target_xy]), cells)
  print(target_cell)


  print('ncells', len(cells))
  target_cell[0].refine(nsub=[4, 4])
  cells = root.get_leaves()
  print('ncells', len(cells))
  ncells = len(cells)
  xy = hyper.HCell.get_centroids(cells)
  F = branin.evaluate(xy[:, 0], xy[:, 1])
  rb = RBFInterpolator(xy, F, kernel="quintic")
  F_est = rb(xy_test)
  error = F_est - F_test
  error = np.absolute(error)
  print('err', error[idx[0]], error[idx[-1]])
  print(np.max(error))


  idx = np.argsort(error)
  #print(error[idx[0]], error[idx[-1]])

def t2():

  np.random.seed(0)
  ntest = 100000
  xy_test = np.random.random((ntest, 2))
  xy_test[:, 0] *= 15.0
  xy_test[:, 0] += -5.0
  xy_test[:, 1] *= 15.0
  F_test = branin.evaluate(xy_test[:, 0], xy_test[:, 1])


  root = hyper.HCell([-5.0, 0.0], end=[10.0, 15.0])

  root.refine(nsub=[12, 12])
  cells = root.get_leaves()

  ncells = len(cells)
  #xy = np.zeros((ncells, 2))
  xy = hyper.HCell.get_centroids(cells)

  F = branin.evaluate(xy[:, 0], xy[:, 1])

  rb = RBFInterpolator(xy, F, kernel="quintic")

  #rb(np.array([[1.0, 1.0]]))



  for iter in range(100):
    print('iter', iter)
    F_est = rb(xy_test)
    error = F_est - F_test
    error = np.absolute(error)
    print('  max err', np.max(error))

    idx = np.argsort(error)
    print('  min/max err', error[idx[0]], error[idx[-1]], 'npoints', xy.shape[0])
    print('  xy_test: coor/idx', xy_test[ idx[-1], : ], idx[-1])



    ns = len(idx)
    idx_ = set( idx[ns-1:] ) # flattern as set to remove dups
    print('n points refine', len(idx_))
    """
    target_xy = xy_test[ idx[-1], : ]
    target_cell = point_in_cell(np.array([target_xy]), cells)
    print('  ', target_cell)


    print('  ncells (before refine)', len(cells))
    target_cell[0].refine(nsub=[4, 4])
    cells = root.get_leaves()
    print('  ncells (after refine)', len(cells))
    ncells = len(cells)
    """

    target_xy = xy_test[ list(idx_), : ] # use list(idx_) to enable slicing
    target_cell = point_in_cell( target_xy, cells)
    print('n cells to refine', len(target_cell))
    target_cell = set(target_cell)
    print('n unique cells refine', len(target_cell))

    print('  ncells (before refine)', len(cells))
    """
    target_cell[0].refine(nsub=[4, 4])
    """
    for tc in target_cell:
        tc.refine(nsub=[2, 2])

    cells = root.get_leaves()
    print('  ncells (after refine)', len(cells))
    ncells = len(cells)


    xy = hyper.HCell.get_centroids(cells)
    F = branin.evaluate(xy[:, 0], xy[:, 1])
    rb = RBFInterpolator(xy, F, kernel="quintic")


  xy = hyper.HCell.get_centroids(cells)
  np.savetxt('xy.gp',xy)


if __name__=="__main__":
  np.random.seed(10)
  npoints = 2314
  xy = np.random.random((npoints, 2))
  xy[:, 0] *= 15.0
  xy[:, 0] += -5.0
  xy[:, 1] *= 15.0
  F = branin.evaluate(xy[:, 0], xy[:, 1])
  rb = RBFInterpolator(xy, F, kernel="quintic")

  np.random.seed(0)
  ntest = 100000
  xy_test = np.random.random((ntest, 2))
  xy_test[:, 0] *= 15.0
  xy_test[:, 0] += -5.0
  xy_test[:, 1] *= 15.0
  F_test = branin.evaluate(xy_test[:, 0], xy_test[:, 1])

  F_est = rb(xy_test)
  error = F_est - F_test
  error = np.absolute(error)
  print('npoints', npoints, 'max err', np.max(error))


  t2()
