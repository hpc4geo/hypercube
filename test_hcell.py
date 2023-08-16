
from hcell import *

# Build a 1-D domain consisting of a single HCell.
def test1():
  root = HCell([0.0], width=[1.0])
  print('bb', root.get_bounding_box())
  print('centroid', root.get_centroid())

  root.refine(nsub=[2])
  print('root', root)
  print('leaves', root.children)

  print('--')
  cells = root.get_leaves()
  for c in cells:
    print(c)
    print('bb', c.get_bounding_box()[0])
    print('centroid', c.get_centroid())

  """
  print('--')
  cells = root.get_cells()
  for c in cells:
    print(c)
  """


# Build a 2-D domain consisting of a single HCell.
def test2():
  root = HCell([0.0, 0.0], width=[1.0, 1.0])
  print('bb', root.get_bounding_box())
  print('centroid', root.get_centroid())

  print('face centroid[x]', root.get_face_centroids(0))
  print('face centroid[y]', root.get_face_centroids(1))

  print('vertices', root.get_vertices())

  cc = root.get_leaves()
  print('nleaves', len(cc)) # should be 0
  cc = root.get_cells()
  print('ncells', len(cc)) # should be 1

  root.refine(nsub=[4]*2)

  print('--')
  cells = root.get_leaves()
  for c in cells:
    print('bb', c.get_bounding_box(), 'centroid', c.get_centroid())

  cc = root.get_leaves()
  print('nleaves', len(cc)) # should be 4 * 4 = 16
  cc = root.get_cells()
  print('ncells', len(cc)) # should be 4 * 4 + 1 = 17


# Build a 2-D domain consisting of 4x3 HCell objects.
def test_mesh():
  cells = HCell.build_base_mesh([0.0, 0.0], [1.0, 2.0], [4,3])
  for c in cells:
    print('bb', c.get_bounding_box(), 'centroid', c.get_centroid())
  print(HCell.get_centroids(cells))

  leaves = HCell.collect_leaves(cells)
  print('leaves', leaves)

  cc = HCell.collect_leaves(cells)
  print('nleaves', len(cc)) # should be 0
  cc = HCell.collect_cells(cells)
  print('ncells', len(cc)) # should be 4 * 3 = 12

  cells[7].refine(nsub=[2,2])

  cc = HCell.collect_leaves(cells)
  print('nleaves', len(cc)) # should be 2 * 2 = 4
  cc = HCell.collect_cells(cells)
  print('ncells', len(cc)) # should be 4 * 3 + 2 * 2 = 16


# Build a 11-D domain consisting of 4 HCell objects in each i dimension.
def test_deep(dim=11):
  # Time starts to get long (> 30 sec) if dim > 11

  import time as time

  logs = list()

  t0 = time.perf_counter()
  cells = HCell.build_base_mesh([0.0]*dim, [1.0]*dim, [4]*dim)
  logs.append(['build_base_mesh', time.perf_counter() - t0])
  print('ncells', len(cells))

  # arbitrarily refine cell 10
  t0 = time.perf_counter()
  cells[10].refine(nsub=[2]*dim)
  logs.append(['refine-1', time.perf_counter() - t0])

  # arbitrarily refine leaf 20
  t0 = time.perf_counter()
  leaves = HCell.collect_leaves(cells)
  logs.append(['collect_leaves-1', time.perf_counter() - t0])
  print('nleaves', len(leaves))

  t0 = time.perf_counter()
  leaves[20].refine(nsub=[2]*dim)
  logs.append(['refine-2', time.perf_counter() - t0])
  t0 = time.perf_counter()
  leaves = HCell.collect_leaves(cells)
  logs.append(['collect_leaves-2', time.perf_counter() - t0])
  print('nleaves', len(leaves))

  # Dump cells and load from file
  t0 = time.perf_counter()
  HCell.dump(cells)
  logs.append(['dump', time.perf_counter() - t0])

  t0 = time.perf_counter()
  cells_ = HCell.load()
  logs.append(['load', time.perf_counter() - t0])

  leaves = HCell.collect_leaves(cells_)
  print('ncells', len(cells_))
  print('nleaves', len(leaves))

  # Flatten, dump, load cells
  cells_ = HCell.collect_leaves(cells)
  HCell.reset_all(cells_)
  print('ncells', len(cells_))
  HCell.dump(cells_)
  cells_ = HCell.load()
  leaves = HCell.collect_leaves(cells_)
  print('ncells', len(cells_))
  print('nleaves', len(leaves))

  print(logs)

if __name__=="__main__":
  #test1()
  #test2()
  test_mesh()
  #test_deep(dim=12)
