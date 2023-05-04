import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation



def complete_disp(bc_array, nodes, sol, ndof_node=2):
    """
    Fill the displacement vectors with imposed and computed values.

    bc_array : ndarray (int)  [node_idx][0:2]
        Indicates if the nodes has any type of boundary conditions
        applied to it.
    sol : ndarray (float)
        Array with the computed displacements.
    nodes : ndarray (float)
        Array with number and nodes coordinates
    ndof_node : int
        Number of degrees of freedom per node.

    Returns
    -------
    sol_complete : (nnodes, ndof_node) ndarray (float)
      Array with the displacements.
    """
    num_node= bc_array.shape[0]
    node_dofs= bc_array.shape[1]
    sol_complete = np.zeros((num_node,node_dofs))

    for i in range(num_node): # loop node 
        for j in range(node_dofs): # loop node dofs
            dof_g = int(bc_array[i][j])
            if dof_g != -1:
                sol_complete[i][j] = sol[dof_g]
            
    return sol_complete

def mesh2tri(nodes, elements):
    """Generate a  matplotlib.tri.Triangulation object from the mesh

    Parameters
    ----------
    nodes : ndarray (float)
      Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
      Array with the node number for the nodes that correspond to each
      element.

    Returns
    -------
    tri : Triangulation , each tri indicate by element(i,j,k)
        An unstructured triangular grid consisting of npoints points
        and ntri triangles.

    """
    coord_x = nodes[:,1]
    coord_y = nodes[:,2]
    tris = []
    for elem in elements:
        eleType = elem[1]
        if eleType == 1:
            tris.append(elem[[3,4,5]])
            tris.append(elem[[5,6,3]])
        if eleType == 2:
            pass
    return Triangulation(coord_x,coord_y,tris)



def plot_node_field(field, nodes, elements, plt_type="contourf", levels=12,
                    savefigs=False, title=None, figtitle=None,
                    filename=None):
    """Plot the nodal displacement using a triangulation

    Parameters
    ----------
    field : ndarray (float)
          Array with the field to be plotted. The number of columns
          determine the number of plots.
    nodes : ndarray (float)
        Array with number and nodes coordinates
        `number coordX coordY BCX BCY`
    elements : ndarray (int)
        Array with the node number for the nodes that correspond
        to each  element.
    plt_type : string (optional)
        Plot the field as one of the options: ``pcolor`` or
        ``contourf``.
    levels : int (optional)
        Number of levels to be used in ``contourf``.
    savefigs : bool (optional)
        Allow to save the figure.
    title : Tuple of strings (optional)
        Titles of the plots. If not provided the plots will not have
        a title.
    figtitle : Tuple of strings (optional)
        Titles of the plotting windows. If not provided the
        windows will not have a title.
    filename : Tuple of strings (optional)
        Filenames to save the figures. Only used when `savefigs=True`.
        If not provided the name of the figures would be "outputk.pdf",
        where `k` is the number of the column.
    """
    tri = mesh2tri(nodes, elements)

    if len(field.shape) == 1:
        nfields = 1
    else:
        _, nfields = field.shape
    if title is None:
        title = ["" for cont in range(nfields)]
    if figtitle is None:
        figs = plt.get_fignums()
        nfigs = len(figs)
        figtitle = [cont + 1 for cont in range(nfigs, nfigs + nfields)]
    if filename is None:
        filename = ["output{}.pdf".format(cont) for cont in range(nfields)]
    for cont in range(nfields):
        if nfields == 1:
            current_field = field
        else:
            current_field = field[:, cont]
        plt.figure(figtitle[cont])

        plt.triplot(tri, 'ro-', lw = 1)
        plt.show()
        
        #tri_plot(tri, current_field, title=title[cont], levels=levels,
        #         plt_type=plt_type, savefigs=savefigs,
        #         filename=filename[cont])
        if savefigs:
            plt.savefig(filename[cont])


def tri_plot(tri, field, title="", levels=12, savefigs=False,
             plt_type="contourf", filename="solution_plot.pdf"):
    """Plot contours over triangulation

    Parameters
    ----------
    tri : ndarray (float)
        Array with number and nodes coordinates:
        `number coordX coordY BCX BCY`
    field : ndarray (float)
        Array with data to be plotted for each node.
    title : string (optional)
        Title of the plot.
    levels : int (optional)
        Number of levels to be used in ``contourf``.
    savefigs : bool (optional)
        Allow to save the figure.
    plt_type : string (optional)
        Plot the field as one of the options: ``pcolor`` or
        ``contourf``
    filename : string (optional)
        Filename to save the figures.
    """
    if plt_type == "pcolor":
        disp_plot = plt.tripcolor
    elif plt_type == "contourf":
        disp_plot = plt.tricontourf
    disp_plot(tri, field, levels, shading="gouraud")
    plt.title(title)
    plt.colorbar(orientation='vertical')
    plt.axis("image")
    plt.show()
    if savefigs:
        plt.savefig(filename)

