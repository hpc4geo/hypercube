

import numpy as np
import itertools as itools

class HCell:

  def __init__(self, start, width=None, end=None, nsub=None, parent=None):
    """
    Input
    -----
      - start: iterable
        Minimum coordinate value in each dimension.
      - width: iterable
        Length of cell in each dimension.
      - end: iterable (optional)
        Maximum coordinate value in each dimension.
      - nsub: iterable (optional)
        Sub-division factor in each dimension.
      - parent: HCell (optional)
        Reference to a cell that defines self as a child cell (i.e. `self in parent->children == True`).
    Exception
    ----------
      - ValueError: Occurs if the dimension of `start` does not match dimension of `end`.
      - ValueError: Occurs if the cell width in any dimension is <= 0.
    """
    self.ndim = len(start)
    self.corner = np.array(start)
    if width is None:
      self.width = np.array(end) - self.corner
    else:
      self.width = np.array(width)
    if np.min(self.width) <= 0.0:
      raise ValueError("Cell width in each directio must be positive")
    if self.corner.shape[0] != self.width.shape[0]:
      raise ValueError("Dimension of `corner` and `width` must be equal")

    self.nsub = None
    self.parent = None
    self.children = None

    if nsub is not None: self.nsub = nsub
    if parent is not None: self.parent = parent

  def get_centroid(self):
    """
    Compute the n-dim centroid of the cell.

    Output
    ------
      - xc: numpy.ndarray, shape = (ndim)
    """
    xc = np.array(self.corner) # force a copy of corner
    xc += 0.5 * self.width
    return xc


  def get_bounding_box(self):
    """
    Compute the n-dim bounding box of the cell.

    Output
    ------
      - x0: numpy.ndarray, shape = (ndim)
        Minimum coordinate values.
      - x1: numpy.ndarray, shape = (ndim)
        Maximum coordinate values.
    """
    x0 = np.array(self.corner) # force a copy of corner
    x1 = x0 + self.width
    return x0, x1


  def get_face_centroids(self, dir):
    """
    Compute the n-dim face centroids of the cell.

    Input
    -----
      - dir: int
        Defines the faces along which centroids will be computed.
    Output
    ------
      - fc: numpy.ndarray, shape = (2, ndim)
        Face centroid associated with the coordinates `self.get_bounding_box()[0][dir]`.
        Face centroid associated with the coordinate `self.get_bounding_box()[1][dir]`.
    Exception
    ---------
      - ValueError: Occurs if `dir` >= self.ndim
    """
    if dir < 0 or dir >= self.ndim:
      raise ValueError("`dir` must in the range [0," + str(self.ndim-1)+"]. Value provided was " + str(dir))
    xmin, xmax = self.get_bounding_box()
    centroid = self.get_centroid()
    fc = np.zeros((2, self.ndim))
    c0 = np.array(centroid)
    c0[dir] = xmin[dir]
    c1 = np.array(centroid)
    c1[dir] = xmax[dir]
    fc[0, :] = c0[:] # Min value in direction `dir`
    fc[1, :] = c1[:] # Max value in direction `dir`
    return fc

  def _descend(self, cells):
    if self.children is None:
      if self.parent is not None:
        cells.append(self)
        return cells
    else:
      for c in self.children:
        c._descend(cells)
    return cells


  def get_vertices(self):
    """
    Compute the vertices of the cell.

    Output
    ------
      - vcoor: numpy.ndarray, shape = (2**ndim, ndim)
        All vertices associated with the cell.
    """
    bbmin, bbmax = self.get_bounding_box()
    coor = list()
    for k in range(self.ndim):
      x1d = [bbmin[k], bbmax[k]]
      coor.append(x1d)
    coords_prod = itools.product(*coor)

    vcoor = np.zeros((2**self.ndim, self.ndim))
    k = 0
    for c in coords_prod:
      vcoor[k, :] = c[:]
      k += 1
    return vcoor


  def get_leaves(self):
    """
    Collects all leaves into a list.
    Leaves are defined as cells which do not have any children.

    Output
    ------
      - cells: list
        All descendents of self which do not have children.
    """
    cells = list()
    self._descend(cells)
    return cells


  def _descend_all(self, cells):
    if self.children is None:
      return cells
    else:
      for c in self.children:
        cells.append(c)
        c._descend_all(cells)
    return cells

  def get_cells(self):
    """
    Collects all cells into a list.

    Output
    ------
      - cells: list
        Contains self and all descendants of self. `cells` will contain unique cells.
    """
    cells = list()
    cells.append(self)
    self._descend_all(cells)
    return cells


  def coarsen(self):
    """
    Coarsen a cell. This simply removes any child cells.

    Exception
    ---------
      - RuntimeError: Occurs if child cells do not exist.
    """
    if self.children is None:
      raise RuntimeError('Cannot coarsen a cell without children')
    self.children = None


  def refine(self, nsub=None):
    """
    Refine a cell. Refinement is defined by sub-dividing a cell in each dimension.

    Input
    -----
      - nsub: iterable (optional)
        The number of sub-divisions in each dimension a cell should be split into.
    Output
    ------
      - child_cells: list
        The new cells created due to refinement.
    Exception
    ---------
      - RuntimeError: Occurs if child cells are present.
      - ValueError: Occurs if nsub is None.
    """
    if self.children is not None:
      raise RuntimeError('Cannot refine a cell previously refined')

    if nsub is None:
      nsub = self.nsub
    if nsub is None:
      raise ValueError('Must provide value for `nsub` for root cells')

    w = np.array(self.width)
    sizes = np.array(nsub)
    w /= sizes
    #print('[refine] w', w)

    #print('coord')
    coords = list()
    for k in range(self.ndim):
      x1d = np.linspace(0.0, self.width[k], nsub[k], endpoint=False)
      x1d += self.corner[k]
      #print('[refine] p', k, 'c0', x1d)
      coords.append(x1d)
    #print('/coord')

    #print('prod')
    coords_prod = itools.product(*coords)
    #print('/prod')

    """
    print('child')
    child_cells = list()
    for c0 in coords_prod:
      #print('[refine] child', c0)
      c = HCell(c0, width=list(w))
      c.nsub = nsub
      c.parent = self
      child_cells.append(c)
    print('/child',len(child_cells))
    """

    #print('child')
    # list comprehension
    w_ = list(w)
    child_cells = [ HCell(c0, width=w_, nsub=nsub, parent=self) for c0 in coords_prod ]
    #print('/child',len(child_cells))

    """
    #print('child')
    # generator
    w_ = list(w)
    child_cells = (HCell(c0, width=w_, nsub=nsub, parent=self) for c0 in coords_prod)
    #for c in child_cells:
    #  c.nsub = nsub
    #  c.parent = self
    #print('/child',len(child_cells))
    """

    self.children = child_cells
    return child_cells


  def reset(self):
    """
    Remove inheritied information associated with refinement.
    """
    self.nsub = None
    self.parent = None
    self.children = None


  @classmethod
  def get_centroids(cls, cells):
    """
    Compute cell centroids for a collection of cells.

    Input
    -----
      - cells: HCell or iterable
        A single, or a collection of, HCell objects.
    Output
    ------
      - coor: numpy.ndarray, shape = (ncells, ndim)
        Centroid coordinate of all cells in `cells`.
    """
    if isinstance(cells, HCell):
      cells = list(cells)
    ndim = cells[0].ndim
    ncells = len(cells)
    coor = np.zeros((ncells, ndim))
    for i in range(ncells):
      coor[i, :] = cells[i].get_centroid()
    return coor


  @classmethod
  def reset_all(cls, cells):
    """
    Applies `reset()` to a collection of cells.

    Input
    -----
      - cells: HCell or iterable
        A single, or a collection of, HCell objects.
    """
    if isinstance(cells, HCell):
      cells = list(cells)
    for c in cells:
      c.nsub = None
      c.parent = None
      c.children = None


  @classmethod
  def collect_cells(cls, cells):
    """
    Gathers all cells and their descendents from a collection of HCell objects.

    Input
    -----
      - cells: HCell or iterable
        A single, or a collection of, HCell objects.
    Output
    ------
      - all: list
        Collection of all cells. This includes parents and children.
    """
    if isinstance(cells, HCell):
      cells = list(cells)
    all = list()
    for c in cells:
      all_ = c.get_cells()
      all += all_
    return all


  @classmethod
  def collect_leaves(cls, cells):
    """
    Gathers all leaf cells from a collection of HCell objects.

    Input
    -----
      - cells: HCell or iterable
        A single, or a collection of, HCell objects.
    Output
    ------
      - all: list
        Collection of leaf cells.
    """
    if isinstance(cells, HCell):
      cells = list(cells)
    all = list()
    for c in cells:
      all_ = c.get_leaves()
      all += all_
    return all


  @classmethod
  def build_base_mesh(cls, start, end, nsub):
    """
    Creates a domain sub-divided into M_i pieces in each i dimension.

    Input
    -----
      - start: iterable
        Minimum coordinate value in each dimension of the domain.
      - end: iterable
        Maximum coordinate value in each dimension of the domain.
      - nsub: iterable
        Number of cells in each dimension.
    Output
    ------
      - cells: list
        Collection of root cells.
    """
    root = HCell(start, end=end)
    #print('refine')
    root.refine(nsub=nsub)
    #print('/refine')
    #print('leaves')
    cells = root.get_leaves()
    #print('/leaves')
    HCell.reset_all(cells)
    return cells


  @classmethod
  def dump(cls, cells, fname='hcell.pkl'):
    """
    Write a single or a collection of HCell objects to file.

    Input
    -----
      - cells: HCell or iterable
        A single, or a collection of, HCell objects.
      - fname: str (optional)
        Filename of the pickle file.
    """
    import pickle as pkl

    with open(fname, 'wb') as fp:
      pkl.dump(cells, fp)

  @classmethod
  def load(cls, fname='hcell.pkl'):
    """
    Load a single or a collection of HCell objects from file.

    Input
    -----
      - cells: HCell or iterable
        A single, or a collection of, HCell objects.
      - fname: str (optional)
        Filename of the pickle file.
    Output
    ------
      - cells: list
        Collection of HCell objects.
    """
    import pickle as pkl

    cells = None
    with open(fname, 'rb') as fp:
      cells = pkl.load(fp)
    return cells
